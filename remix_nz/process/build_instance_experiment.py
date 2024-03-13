# %% 
from will_build_instance import *
from remix.framework.tools.gdx import GDXEval
from remix.framework.tools.plots import plot_network, plot_choropleth
import geopandas as gpd
import warnings

# %% [markdown]
# Define useful lists

yrs2run=[2020,2035,2050] # years to be optimised [2020,2030,2040,2050] #
h2_annualdemand=0 # if it is one it takes the annual demand
yrs_demand= yrs2run
yrs_mentioned= [2000, 2020, 2025, 2030, 2035, 2040, 2045, 2050] # must include all years that data is provided for in the model
# Demand files available for different scenarios
files_lst=["01-battery-distributed","01-h2-distributed", "02-battery-overnight", "02-h2-overnight", "03-battery-recharging", "03-h2-recharging", "04-battery-solar", "04-h2-solar"]
files_lst=["02-battery-overnight"]
h2_lst=["01-h2-distributed", "02-h2-overnight",  "03-h2-recharging",  "04-h2-solar"]
elec_lst=["01-battery-distributed", "02-h2-overnight", "03-battery-recharging", "04-battery-solar"]
files_lst=elec_lst
# OJO LO ESTOY CORRIENDO SIN H2

# Path definition
path_base = "C:/Local/REMix" #"C:/Local/REMix/projects/remix_nz" 
path_input = f"{path_base}/remix_nz/input"
path_output = f"{path_base}/remix_nz/output/will" 
path_demand = f"{path_input}/demand/will"       # eletricity, h2
path_profiles = f"{path_input}/profiles"        # renewables, hydro
path_geo = f"{path_input}/shapefiles"           # geojson
path_brownfield = f"{path_input}/brownfield" 

#path_figures = f"{path_output}/figures"        # not in use yet

# %% [markdown]
# Build data

def build_from_list(lst=files_lst):
    nodes_lst=["NIS","AKL","WTO","TRN","BOP","HBY","CEN","WEL","NEL","CAN","OTG" ]
    s3 = time.perf_counter()
    for index, item in enumerate(lst): #enumerate provides both the index and the item in each iteration
        print(f"Building {index+1} of {len(lst)}: {item}")
        yrs_str='-'.join([str(item) for item in yrs2run])
        demand_file=lst[index] 
        yrs_str='-'.join([str(item) for item in yrs2run])
        case_name=f"{demand_file}_{yrs_str}"
        geofile="11regionsNZ.geojson"

        # output
        data_dir = Path(f"{path_output}/{case_name}/data").mkdir(parents=True, exist_ok=True)
        results_dir = Path(f"{path_output}/{case_name}/result").mkdir(parents=True, exist_ok=True)
        print(f"----------Creating data for {case_name}")

        if __name__ == "__main__":
            # Create instance
            s1 = time.perf_counter()
            m = Instance()

            add_nodes(m)
            add_demand(m)
            add_scope(m)

            # renewables
            add_renewables(m)
            add_geothermal(m)
            add_hydro(m)

            # storage
            add_lithium_batteries(m)
            
            # conventional
            add_thermal(m)
            add_gas_turbines(m)

            # hydrogen
            if files_lst==h2_lst:
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
    # Create data
    e3 = time.perf_counter()
    d3=time.strftime("%Hh %Mm %Ss", time.gmtime(e3-s3))
    print(f"Writing dataset for {len(lst)} files took {d3}.")
#build_from_list()

# %% [markdown]
# Run the model

def run_list(lst=files_lst):
    s2 = time.perf_counter()
    for index, item in enumerate(lst): #enumerate provides both the index and the item in each iteration
        print(f"Run {index+1} of {len(lst)}: {item}")
        yrs_str='-'.join([str(item) for item in yrs2run])
        demand_file=files_lst[index] 
        case_name=f"{demand_file}_{yrs_str}"
        output_dir = Path(f"{path_output}/{case_name}")
        data_dir = Path(f"{output_dir}/data")
        result_dir = Path(f"{output_dir}/result")
        if not data_dir.exists():
            raise IOError("You need to build the data!")
        m = Instance.from_path(data_dir)
        # running GAMS from Python script
        s1 = time.perf_counter()
        m.run(
            resultdir = result_dir,
            resultfile=case_name,
            lo=4,
            timeres=1,
            names=1,
            roundts=1,
            iis=1,
        )
        print(os. getcwd())
        # end counter 
        e1 = time.perf_counter()
        d1=time.strftime("%Hh %Mm %Ss", time.gmtime(e1-s1))
        print(f"------------- Running {item} took {d1}.")
    # end counter 
    e2 = time.perf_counter()
    d2=time.strftime("%Hh %Mm %Ss", time.gmtime(e2-s2))
    print(f"------------- Running all cases took {d2}.")
#run_list()

# %% [markdown]
# Build data
def build_data(demand_file,h2="no-h2"):
    nodes_lst=["NIS","AKL","WTO","TRN","BOP","HBY","CEN","WEL","NEL","CAN","OTG" ]
    yrs_str='-'.join([str(item) for item in yrs2run])
    case_name=f"{demand_file}_{yrs_str}"
    geofile="11regionsNZ.geojson"

    # output
    data_dir = Path(f"{path_output}/{case_name}/data").mkdir(parents=True, exist_ok=True)
    results_dir = Path(f"{path_output}/{case_name}/result").mkdir(parents=True, exist_ok=True)
    print(f"----------Creating data for {case_name}")

    if __name__ == "__main__":
        # Create instance
        s1 = time.perf_counter()
        m = Instance()

        add_nodes(m)
        add_demand(m)
        add_scope(m)

        # renewables
        add_renewables(m)
        add_geothermal(m)
        add_hydro(m)

        # batteries
        add_lithium_batteries(m)
        
        # conventional
        add_thermal(m)
        add_gas_turbines(m)

        if h2=="yes-h2":
            # hydrogen
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

# %% [markdown]
# Run the model
def case_run(demand_file):
    yrs_str='-'.join([str(item) for item in yrs2run])
    case_name=f"{demand_file}_{yrs_str}"
    output_dir = Path(f"{path_output}/{case_name}")
    data_dir = Path(f"{output_dir}/data")
    result_dir = Path(f"{output_dir}/result")
    if not data_dir.exists():
        raise IOError("You need to build the data!")
    m = Instance.from_path(data_dir)
    s1 = time.perf_counter()
    # running GAMS from Python script
    m.run(
        resultdir = result_dir,
        resultfile=case_name,
        lo=4,
        timeres=1,
        names=1,
        roundts=1,
        iis=1,
    )
    print(os. getcwd())
    e1 = time.perf_counter()
    d1=time.strftime("%Hh %Mm %Ss", time.gmtime(e1-s1))
    print(f"------------- Running {demand_file} took {d1}.")

# %% [markdown]
# Export results to non GAMS/python users
def results_to_csv(demand_file):
    case_name=f"{demand_file}_{yrs_str}"
    # read in the output `*.gdx` file from the optimization in GAMS
    path_result=Path(f"{path_output}/{case_name}/result")
    path_csv=Path(f"{path_result}/csv_{case_name}")#.mkdir(parents=True, exist_ok=True)
    results = GDXEval(f"{path_result}/{case_name}.gdx")

    # Converter capacities
    converter_caps = results["converter_caps"]   # convert converter capacities to a Pandas DataFrame
    converter_caps = converter_caps[converter_caps > 0.01].dropna()  # remove all capacities with less than 10 MW
    converter_caps.to_csv(f"{path_csv}/{case_name}_converter_caps.csv")
    # Example of slicing:
    # accNodesModel=all, accYears=2050 only, techs=all, commodit=only Elec, capType=total
    # caps.loc[idx[:, "2050", :, "Elec", "total"], :].round(2) 

    # Storage capacities
    storage_caps = results["storage_caps"]   # convert storage capacities to a Pandas DataFrame
    storage_caps = storage_caps[storage_caps > 0.01].dropna()  # remove all capacities with less than 10 MW
    storage_caps.to_csv(f"{path_csv}/{case_name}_storage_caps.csv")

    # Commodity balance (8760 per year) too larg for excel
    commodity_balance = results["commodity_balance"]   # convert commodity_balance to a Pandas DataFrame
    commodity_balance.to_csv(f"{path_csv}/{case_name}_commodity_balance.csv")

    # Commodity balance annual
    commodity_balance_annual = results["commodity_balance_annual"]   # convert commodity_balance to a Pandas DataFrame
    commodity_balance_annual.to_csv(f"{path_csv}/{case_name}_commodity_balance_annual.csv")

    # Marginals
    marginals_sourcesink_profile = results["marginals_sourcesink_profile"]   # convert marginals to a Pandas DataFrame
    marginals_sourcesink_profile.to_csv(f"{path_csv}/{case_name}_marginals_sourcesink_profile.csv")

    # Accounting indicators
    indicator_accounting = results["indicator_accounting"]   # convert  to a Pandas DataFrame
    indicator_accounting.to_csv(f"{path_csv}/{case_name}_indicator_accounting.csv")



for file in files_lst:
    s1 = int(time.perf_counter())
    build_data(file)
    e1 = time.perf_counter()
    d1 = time.strftime("%Hh %Mm %Ss", time.gmtime(e1-s1))
    print(f"Driting dataset took {d1}.")

    s2 = int(time.perf_counter())
    case_run(file)
    e2 = time.perf_counter()
    d2 = time.strftime("%Hh %Mm %Ss", time.gmtime(e2-s2))
    print(f"Running the model took {d2}.")

    #s3 = int(time.perf_counter())
    #results_to_csv(file)
    e3 = time.perf_counter()
    #d3 = time.strftime("%Hh %Mm %Ss", time.gmtime(e3-s3))
    d4 = time.strftime("%Hh %Mm %Ss", time.gmtime(e3-s1))

    print(f"Driting dataset took {d1}.")
    print(f"Running the model took {d2}.")
    #print(f"Exporting to csv files took {d3}.")
    print(f"----------- Total time for {file}: {d4}.")

    
    
# %%
