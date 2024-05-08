
# %% [markdown]
# ## Part A: setting up the model

# Import dependencies
import os
import time
import json
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from remix.framework.api.instance import Instance
from remix.framework.tools.plots import plot_network, plot_choropleth
import networkx as nx
import gams.transfer as gt
from IPython.display import display
idx = pd.IndexSlice

# %% [markdown]
# ### Define the years to run the optimisation and the demand file
# The demand file and the years run determine the name of the case and its results
will_elec = ["00-test-elec","01-battery-distributed", "02-battery-overnight", "03-battery-recharging", "04-battery-solar"]
will_h2 = ["01-h2-distributed", "02-h2-overnight", "03-h2-recharging", "04-h2-solar"]
sdewes_ap = ["base", "high", "med", "ev-med", "low"]
mbie=["base","h2"]
    
scenario_dict = {       
    "will": [will_h2, [2020, 2035, 2050]],
    "sdewes-ap": [sdewes_ap, [2020, 2030, 2040, 2050]],
    "mbie": [mbie, [2020, 2030, 2040, 2050]]
}
group_name="will"
files_lst = scenario_dict[group_name][0]
yrs_sel = scenario_dict[group_name][1] # [2020, 2025, 2030, 2035, 2040, 2045, 2050]
yrs_str='-'.join([str(item) for item in yrs_sel])
yrs_to_calc = [2020, 2025, 2030, 2035, 2040, 2045, 2050]
indx=3

# Define paths/directories
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


def add_scope(m):
    df = pd.DataFrame(
        [
            ["AKL", "AKL", 1],  
            ["BOP", "BOP", 1],  
            ["NEL", "NEL", 1],  
            ["NIS", "NIS", 1],  
            ["OTG", "OTG", 1], 
            ["TRN", "TRN", 1],  
            ["WEL", "WEL", 1], 
            ["WTO", "WTO", 1],  
            ["CAN", "CAN", 1],
            ["CEN", "CEN", 1],  
            ["HBY", "HBY", 1],  
        ]
    )
    df.columns = ["nodesData", "nodesModel", "aggregate"]
    df.set_index(["nodesData", "nodesModel"], inplace=True)
    df["aggregate"] = ""
    df.columns = [""]

    m.map.add(df, "aggregatenodesmodel")

    # Get the data and model regions based on the mapping
    # "set_nodesData.dat"
    m.set.add(list(sorted(set(m.map.aggregatenodesmodel.index.get_level_values(0)))), "nodesdata")
    # "set_nodesModel" & "set_nodesModelSel"
    m.set.add(list(sorted(set(m.map.aggregatenodesmodel.index.get_level_values(1)))), "nodesmodel")

    # Set the years to be considered in the model and the years to be optimized
    # "set_years"
    m.set.add(yrs_to_calc, "years")  # must include all years that data is provided for in the model
    # "set_yearsSel"
    m.set.add(yrs_sel, "yearssel")  # years to be optimised

def add_demand(m):
    rename_commodity = {"Electricity": "Elec",
                        "Hydrogen": "H2",
                        "H2-feedstock": "H2",
                        "Natural Gas": "CH4",
                        "Gas": "CH4",
                        "Feedstock Gas": "CH4",
                        "Feedstock methanol": "CH3OH",
                        "Renewable Fuels": "REfuel",
                        }
    

    # FFE extremos profiles
    # note: we added the .round(3) part because we were getting errors de
    #ts_ffe = -1 * pd.read_csv("input/nz_profile_11nodes.csv", index_col=[0, 1, 2, 3]).rename(index=rename_commodity)
    ts_ffe = -1 * pd.read_csv(Path(path_demand).joinpath(f"{demand_file}.csv"), index_col=[0, 1, 2, 3]).rename(index=rename_commodity)
    ts_ffe["type"] = "fixed"
    ts_ffe_fixed = ts_ffe.set_index("type", append=True).round(3)
    m.profile.add(ts_ffe_fixed, "sourcesink_profile")

    ts_ffe_cfg = pd.DataFrame(index=ts_ffe.index)
    ts_ffe_cfg["usesFixedProfile"] = 1
    m.parameter.add(ts_ffe_cfg, "sourcesink_config")

    # Slack for electricity
    # for 2020 (minus 31st dec to attain to 8760 hours in leap year): https://www.emi.ea.govt.nz/Forward%20markets/Reports/0NQPKT?DateFrom=20211008&DateTo=20221007&Maturity=SHORT&_rsdr=L364D&_si=v|3

    slack_annual = ts_ffe_cfg.loc[idx[:, :, "Wholesale", "Elec"], idx[:]]
    slack_annual = slack_annual.rename(index={"Wholesale": "Slack"}, columns={"usesFixedProfile": "upper"})
    slack_annual["upper"] = np.inf
    m.parameter.add(slack_annual, "sourcesink_annualsum")

    slack_cfg = slack_annual
    slack_cfg = slack_cfg.rename(columns={"upper": "usesUpperSum"}).replace(np.inf, 1)
    slack_cfg["usesLowerProfile"] = 1
    # display(slack_cfg)
    m.parameter.add(slack_cfg, "sourcesink_config")

    slack_cost = pd.DataFrame(
        index=pd.MultiIndex.from_product([["SlackCost"], ["global"], m.set.years, ["Slack"], ["Elec"]])
    )
    slack_cost["perFlow"] = 10  # EEX Strom Futures Cap 3.000 EUR/MWh -> 3 MEUR/GWh
    # display(slack_cost)
    m.parameter.add(slack_cost, "accounting_sourcesinkflow")


  
    # # Hydrogen
    # # "sourcesink_annualSum"
    h2_nodes =  m.set.nodesdata
    h2_annual = pd.DataFrame(
        index=pd.MultiIndex.from_product([h2_nodes, yrs_to_calc, ["Demand"], ["H2"]])
    )

    # Adding the new column "upper" to the DataFrame and setting values for the year 2050
    h2_annual["upper"] = 0 # Initializing with  0
    h2_annualdemand=0 # if it is one it takes the annual demand
    if h2_annualdemand==1:
        h2_2050 = {
            'NIS': 1038.15623,
            'AKL': 70793.7597,
            'WTO': 11663.81169,
            'BOP': 11700.05017,
            'HBY': 13060.32836,
            'TRN': 62016.13938,
            'CEN': 11842.51657,
            'WEL': 34654.12968,
            'NEL': 9158.885102,
            'CAN': 39314.18837,
            'OTG': 18806.67305
            }
        for node, value in h2_2050.items():
            h2_annual.loc[(node, 2050, "Demand", "H2"), "upper"] = value
    m.parameter.add(h2_annual, "sourcesink_annualsum")

    h2_cfg = pd.DataFrame(
        index=pd.MultiIndex.from_product([h2_nodes, yrs_to_calc, ["Demand"], ["H2"]])
    )
    h2_cfg["usesUpperSum"] = 1
    h2_cfg["usesUpperProfile"] = 1
   
    # h2_cfg["usesLowerProfile"] = 1
    m.parameter.add(h2_cfg, "sourcesink_config")   

    # Derive region and time scope
    m.set.add(list(ts_ffe.index.get_level_values(0)), "nodesdata")
    m.set.add(list(ts_ffe.index.get_level_values(1)), "years")

# renewables
    
def load_feedin_csv(year, aggregate=False, norm=True):
    inst = pd.read_csv(Path(path_profiles).joinpath("region_statistics_2012.csv"), index_col=[0, 1])
    ts = pd.read_csv(Path(path_profiles).joinpath(f"timeseries_{year}.csv"), index_col=[0, 1, 2])

    #C: we adjusted plus half an hour
    # Split year and time index
    ts["year"] = (pd.DatetimeIndex(ts.index.get_level_values(0)) + pd.Timedelta(hours=0.5)).year
    ts["t_model"] = ((pd.DatetimeIndex(ts.index.get_level_values(0)) + pd.Timedelta(hours=0.5)).dayofyear - 1) * 24 + (pd.DatetimeIndex(ts.index.get_level_values(0)) - pd.Timedelta(hours=1)).hour + 1
    ts = ts.reset_index().set_index(["region", "technology", "year", "t_model"]).drop(columns=["t"])

    if aggregate:
        inst = inst.groupby(["technology"]).sum()
        ts = ts.groupby(["technology", "year", "t_model"]).sum()

    if norm:
        ts = ts.unstack(["year", "t_model"])
        ts = ts.loc[idx[inst.index]]

        ts = pd.DataFrame(
            ts.values / pd.DataFrame(inst["installable_per_region"]).values,
            index=ts.index,
            columns=ts.columns,
        )
        ts = ts.stack(["year", "t_model"])

    # TODO: Add solar multiple to endask / or model solar field and heat receiver separately
    # FIXME: Remove dirty hack
    # C: COMMENTED BECAUSE THERE IS NO CSP YET
    #ts.loc[idx[:, ["csp_parabolic_trough"], :, :]] = ts.loc[idx[:, ["csp_parabolic_trough"], :, :]] * 1.6
    #ts.loc[idx[:, ["csp_solar_tower"], :, :]] = ts.loc[idx[:, ["csp_solar_tower"], :, :]] * 2.25
    ts[ts > 1] = 1
    return ts.round(3)

def add_renewables(m):
    
    re_inst_csv = pd.read_csv(Path(path_profiles).joinpath("region_statistics_2012.csv"), index_col=[0, 1])
    re_nodes = [n for n in m.set.nodesdata]

    re_vintage = [1950,2030,2040,2050]
    year_mapping = {
        2012: 2020,
        # 2011: 2025,
        # 2012: 2030,
        #2012: 2035,
        #2014: 2040,
        #2016: 2045,
        #2012: 2050 #this mapping is only used with different weather years
    }
    re_techs = list(set(re_inst_csv.index.get_level_values(1)))
    pv_techs = [i for i in re_techs if i.startswith("pv")]
    csp_techs = []#[i for i in re_techs if i.startswith("csp")]
    wind_techs = [i for i in re_techs if i.startswith("wind")]
    
    re_tech = pd.DataFrame(index=pd.MultiIndex.from_product([re_techs, re_vintage]))
    re_tech.loc[idx[:, :], "activityUpperLimit"] = 0 # so it is overwritten by the availabity from the timeseriesfiles
    re_tech.loc[idx[pv_techs, [1950]], "lifeTime"] = 35  # years
    re_tech.loc[idx[pv_techs, [2030, 2040, 2050]], "lifeTime"] = 40  # years
    re_tech.loc[idx[csp_techs + wind_techs, [1950]], "lifeTime"] = 27  # years
    re_tech.loc[idx[csp_techs + wind_techs, [2030, 2040, 2050]], "lifeTime"] = 30  # years
    m.parameter.add(re_tech, "converter_techparam")
    
    # capacities
    years_build = list(range(1989, 2021)) #brownfield

    re_caps = pd.DataFrame(index=pd.MultiIndex.from_product([re_nodes, years_build, re_techs]))
    re_caps.index.names = ["region", "years", "technology"]

    re_caps.loc[idx["CAN", [2003], "wind_onshore"], "unitsBuild"] = 0.0005 # GW_el
    re_caps.loc[idx["CAN", [2005], "wind_onshore"], "unitsBuild"] = 0.0001 # GW_el

    re_caps.loc[idx["CEN", [1999], "wind_onshore"], "unitsBuild"] = 31.7  /1000 # GW_el	
    re_caps.loc[idx["CEN", [2004], "wind_onshore"], "unitsBuild"] = 127.05 /1000 # GW_el	
    re_caps.loc[idx["CEN", [2007], "wind_onshore"], "unitsBuild"] = 93 /1000 # GW_el
    re_caps.loc[idx["CEN", [2011], "wind_onshore"], "unitsBuild"] = 48.5 /1000 # GW_el
    re_caps.loc[idx["CEN", [2020], "wind_onshore"], "unitsBuild"] = 221.4 /1000 # GW_el	

    re_caps.loc[idx["OTG", [2007], "wind_onshore"], "unitsBuild"] = 58/1000 # GW_el	
    re_caps.loc[idx["OTG", [2009], "wind_onshore"], "unitsBuild"] = 2.25/1000 # GW_el	
    re_caps.loc[idx["OTG", [2010], "wind_onshore"], "unitsBuild"] = 0.45/1000 # GW_el	
    re_caps.loc[idx["OTG", [2011], "wind_onshore"], "unitsBuild"] = 43.65/1000 # GW_el	
    re_caps.loc[idx["OTG", [2015], "wind_onshore"], "unitsBuild"] = 6.8 /1000 # GW_el
    
    re_caps.loc[idx["NEL", [2010], "wind_onshore"], "unitsBuild"] = 0.75/1000 # GW_el
    re_caps.loc[idx["NEL", [2011], "wind_onshore"], "unitsBuild"] = 1/1000 # GW_el
    re_caps.loc[idx["NEL", [2014], "wind_onshore"], "unitsBuild"] = 0.66/1000 # GW_el

    re_caps.loc[idx["TRN", [2020], "wind_onshore"], "unitsBuild"] = 0.1333 # GW_el

    re_caps.loc[idx["WEL", [1993], "wind_onshore"], "unitsBuild"] = 8.45 /1000# GW_el
    re_caps.loc[idx["WEL", [1996], "wind_onshore"], "unitsBuild"] = 8.45 /1000# GW_el
    re_caps.loc[idx["WEL", [2009], "wind_onshore"], "unitsBuild"] = 143 /1000# GW_el
    re_caps.loc[idx["WEL", [2014], "wind_onshore"], "unitsBuild"] = 71.3 /1000# GW_el	
    re_caps.loc[idx["WTO", [2011], "wind_onshore"], "unitsBuild"] = 64.4/1000 # GW_el

    re_upper = pd.DataFrame(index=pd.MultiIndex.from_product([re_nodes, yrs_to_calc, re_techs]))
    re_upper.index.names = ["region", "years", "technology"]
    re_max_inst = pd.DataFrame(re_inst_csv.div(1e3)["installable_per_region"]).rename(columns={"installable_per_region": "unitsUpperLimit"})
    re_upper = re_upper.join(re_max_inst, on=["region", "technology"], how="outer")
    re_upper.loc[idx[:, [2020], :], "noExpansion"] = 1  # boolean
    re_caps_full = pd.concat([re_caps, re_upper], axis=1).dropna(how="all")

    m.parameter.add(re_caps_full, "converter_capacityparam")


    # activity
    # FIXME: (PART 1): it is ok while we only have 1 weather year
    re_feedin_csv = load_feedin_csv(2012).unstack("t_model").swaplevel(1, 2)  # Load data for the year 2012
    re_feedin = pd.concat([re_feedin_csv.rename(index={2012:y}) for y in yrs_to_calc])
    re_feedin["type"] = "upper"
    re_feedin = re_feedin.set_index("type", append=True)
    re_feedin = re_feedin[re_feedin >= 0.01].dropna(how="all").fillna(0)
    re_feedin = re_feedin.iloc[:, 0:8760]
    re_feedin.columns = [f"t{str(t+1).zfill(4)}" for t in range(8760)]
    re_feedin = re_feedin.sort_index(level=["region", "technology", "year"])
    m.profile.add(re_feedin, "converter_activityprofile")

    # coefficients
    re_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [
                re_techs,
                re_vintage,
                ["Powergen", "Heatgen"],
                ["Elec", "Heat_CSP"],
            ]
        )
    )
    re_coef.loc[idx[wind_techs + pv_techs, :, "Powergen", "Elec"], "coefficient"] = 1
    re_coef.loc[idx[csp_techs, :, "Heatgen", "Heat_CSP"], "coefficient"] = 1
    m.parameter.add(re_coef, "converter_coefficient")

    re_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product([["Invest", "OMFix"], ["global"], re_techs, re_vintage])
    ).sort_index()

    # TODO: Update costs for renewable technologies
    # CSP own assumptions based on: https://aip.scitation.org/doi/pdf/10.1063/5.0028883, https://elib.dlr.de/186998/1/SolarPACES_2021_Paper_Dersch_R1.pdf
    #re_acc.loc[idx["Invest", :, "csp_parabolic_trough", :], "perUnitBuild"] = [344.5, 274.7, 230.2, 196.0]  # Child 2019 - Mio EUR per unit
    #re_acc.loc[idx["Invest", :, "csp_solar_tower", :], "perUnitBuild"] = [482, 372, 310, 264]  # Mio EUR per unit

    re_acc.loc[idx["Invest", :, "pv_decentral", :], "perUnitBuild"] = [870, 570, 460, 410]  # DEA2022 PV comm&indust - Mio EUR per unit
    re_acc.loc[idx["Invest", :, "pv_central_fixed", :], "perUnitBuild"] = [560, 380, 320, 290]  # DEA2022 utility scale - Mio EUR per unit
    re_acc.loc[idx["Invest", :, "pv_central_track_azimuth", :], "perUnitBuild"] = [650, 450, 380, 350]  # DEA2022 utility scale (tracking) - Mio EUR per unit

    re_acc.loc[idx["Invest", :, "wind_onshore", :], "perUnitBuild"] = [1330, 1040, 980, 960]  # DEA2022 onshore - Mio EUR per unit
    re_acc.loc[idx["Invest", :, "wind_offshore_foundation", :], "perUnitBuild"] = [2120, 2287, 2168, 2130] # DEA2022 offshore - Mio EUR per unit
    re_acc.loc[idx["Invest", :, "wind_offshore_floating", :], "perUnitBuild"] = 1.2 * np.array([2120, 2287, 2168, 2130])  # DEA2022 offshore + 20% assumption - Mio EUR per unit

    re_acc.loc[idx["Invest", :, pv_techs, [1950]], "amorTime"] = 35  # years
    re_acc.loc[idx["Invest", :, pv_techs, [2030, 2040, 2050]], "amorTime"] = 40  # years
    re_acc.loc[idx["Invest", :, csp_techs + wind_techs, [1950]], "amorTime"] = 27  # years
    re_acc.loc[idx["Invest", :, csp_techs + wind_techs, [2030, 2040, 2050]], "amorTime"] = 30  # years

    re_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1  # binary yes/no
    re_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06  # percent/100
    re_acc.loc[idx["OMFix", :, pv_techs + wind_techs, :], "perUnitTotal"] = (
        re_acc.loc[idx["Invest", :, pv_techs + wind_techs, :], "perUnitBuild"] * 0.02
    )  # Mio EUR per unit
    re_acc.loc[idx["OMFix", :, csp_techs, :], "perUnitTotal"] = (
        re_acc.loc[idx["Invest", :, csp_techs, :], "perUnitBuild"] * 0.015
    )  # Mio EUR per unit

    m.parameter.add(re_acc, "accounting_converterunits")  

def add_geothermal(m):
    geoth_inst_csv = pd.read_csv(Path(path_brownfield).joinpath("power-plant-nz-database.csv"))
    geoth_vintage = [1950]
    geoth_techs = ["geoth"]
    #geoth_nodes = [n for n in m.set.nodesdata if not n.startswith("LNG")]
    geoth_nodes = ["BOP", "NIS", "WTO"]
    geoth_activities = ["Powergen"]

    # techparam 
    geoth_tech = pd.DataFrame(index=pd.MultiIndex.from_product([geoth_techs, geoth_vintage]))
    geoth_tech.loc[idx[:, :], "activityUpperLimit"] = 0.9 # 0 for renewables 
    geoth_tech.loc[idx[:, :], "lifeTime"] = 100  # years, data from: "Financial_Technical assumptions" Ashish 2023 
    m.parameter.add(geoth_tech, "converter_techparam")
    
    df = geoth_inst_csv
    filtered_df = df[df['Type'] == 'Geothermal']
    #geoth_df = filtered_df.groupby(['Node', 'Year_built', 'Techs']).agg({'Capacity_MW': 'sum', 'Avg_Ann_Gen_GWh': 'sum'}).reset_index()            
    grouped_df = filtered_df.groupby(['Node', 'Year_built', 'Techs'])['Capacity_MW'].sum().reset_index()
    grouped_df['Year_built'] = grouped_df['Year_built'].astype(int)

    # capacities
    geoth_cap = (grouped_df
                 .set_index(["Node", "Year_built", "Techs"])
                 .rename(columns={"Capacity_MW": "unitsBuild"})
                 .div(1e3))

    geoth_cap_upper = pd.DataFrame(index=pd.MultiIndex.from_product([geoth_nodes, yrs_to_calc, geoth_techs]))
    geoth_cap_upper.loc[idx[:, :, :], "noExpansion"] = 1  # boolean
    geoth_cap_upper.loc[idx[:, :, :], "unitsUpperLimit"] = 100 

    geoth_cap = pd.concat([geoth_cap, geoth_cap_upper]).groupby(level=[0,1,2]).sum()
    m.parameter.add(geoth_cap, "converter_capacityparam")
    
    # activity data
    geoth_inst = geoth_inst_csv[geoth_inst_csv["Type"] == "Geothermal"]

    geoth_feedin = pd.DataFrame(
        data={
            "region": geoth_inst["Node"],  # Cambiar "Node" por "region"
            "year": np.nan,
            "technology": "Geothermal",  # Agregar "technology" desde el inicio
            "type": "upper",  # Agregar "type" desde el inicio
        },
        columns=["region", "year", "technology", "type"],  # Establecer orden de columnas
    )

    geoth_inst_filtered = geoth_inst[
        (geoth_inst["Avg_Ann_Gen_GWh"] != 0) & (geoth_inst["Capacity_MW"].notnull())
    ]

    # Create the DataFrame efficiently
    data = {}
    for i in range(1, 8761):
        data["t{:04d}".format(i)] = np.nan

    # Calculate plant factors using the filtered DataFrame
    for i in range(1, 8761):
        factor = geoth_inst_filtered["Avg_Ann_Gen_GWh"] / (geoth_inst_filtered["Capacity_MW"]*8760 / 100)
        data["t{:04d}".format(i)] = factor

    geoth_feedin_temp = pd.DataFrame(data=data)
    geoth_feedin_grouped = geoth_feedin_temp.groupby(geoth_inst["Node"]).sum()
    geoth_feedin_grouped.reset_index(inplace=True)
    geoth_feedin_grouped["technology"] = "Geothermal"  
    geoth_feedin_grouped["type"] = "upper"  

    dfs = []
    for year in yrs_to_calc:
        geoth_feedin_copy = geoth_feedin_grouped.copy()
        geoth_feedin_copy["year"] = year
        dfs.append(geoth_feedin_copy)

    # Concatenate the separate dataframes into one
    geoth_feedin = pd.concat(dfs, ignore_index=True)
    geoth_feedin = geoth_feedin.rename(columns={"Node": "region"})
    geoth_feedin = geoth_feedin.set_index(["region", "year", "technology", "type"]).sort_index()  # Establecer índice ahora que las columnas existen
    m.profile.add(geoth_feedin, "converter_activityprofile")
    
    # coefficients    
    geoth_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [geoth_techs, geoth_vintage, geoth_activities, ["Elec"]]
        )
    )

    geoth_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = 1 # GW_el
    m.parameter.add(geoth_coef, "converter_coefficient")

    geoth_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product([["Invest", "OMFix"], ["global"], geoth_techs, geoth_vintage])
    ).sort_index()

    # TODO: Update costs for renewable technologies
    # CSP own assumptions based on: https://aip.scitation.org/doi/pdf/10.1063/5.0028883, https://elib.dlr.de/186998/1/SolarPACES_2021_Paper_Dersch_R1.pdf
    #re_acc.loc[idx["Invest", :, "csp_parabolic_trough", :], "perUnitBuild"] = [344.5, 274.7, 230.2, 196.0]  # Child 2019 - Mio EUR per unit
    #re_acc.loc[idx["Invest", :, "csp_solar_tower", :], "perUnitBuild"] = [482, 372, 310, 264]  # Mio EUR per unit

    geoth_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] = 4970   #, 3610]   # data from: "Financial_Technical assumptions" Ashish 2023   - Mio EUR per unit

    geoth_acc.loc[idx["Invest", :, :, :], "amorTime"] = 100  # years

    geoth_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1  # binary yes/no
    geoth_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06  # percent/100
    geoth_acc.loc[idx["OMFix", :, :, :], "perUnitTotal"] = (
        geoth_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] * 0.08 # data from: "Financial_Technical assumptions" Ashish 2023 
    )  # Mio EUR per unit

    m.parameter.add(geoth_acc, "accounting_converterunits")

def add_hydro(m):

    # "sourcesink_config" (import configuration)
    sourcesink_config = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [m.set.nodesdata, m.set.yearssel, ["Ocean"], ["Water_out"]]
        )
    )
    # fix: negative profile for minimum flow, upper profile of 0
    sourcesink_config.loc[idx[m.set.nodesdata, :, :, :], "usesUpperProfile"] = 1 #we need negative values to get water out of the system
    sourcesink_config.dropna(inplace=True)

    m.parameter.add(sourcesink_config, "sourcesink_config")
    sourcesink_config

    #error from logfile Infeasibility row 'Eq_balance_commodities(tm1,HBY,2020,Water_in)':  0  = -15.39.
    hydro_vintage = [1950]
    hydro_years= [2000]+yrs_to_calc
    hydro_techs = ["Hydro"] # unifying storage (dam) and converter (turbine) in one 
    hydro_nodes = ["BOP", "CAN", "CEN", "HBY", "NEL", "OTG", "WTO"] #[n for n in m.set.nodesdata if not n.startswith("LNG")]
    hydro_activities = ["Power_gen","Spill"] 

	# Converter (turbine)
    conv_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([hydro_techs, hydro_vintage])
    )
    conv_tech.loc[idx[:, :], ["lifeTime"]] = 100
    conv_tech.loc[idx[:, :], ["activityUpperLimit"]] = 1
    m.parameter.add(conv_tech, "converter_techparam")

    conv_cap = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [hydro_nodes, hydro_years, hydro_techs] #fix: maybe years_mentioned
        )
    )
    

    # existing capacities
    conv_cap.loc[idx[["BOP"], [2000], "Hydro"], "unitsBuild"] = 0.17095  # GW_el
    conv_cap.loc[idx[["CAN"], [2000], "Hydro"], "unitsBuild"] = 1.82683  # GW_el
    conv_cap.loc[idx[["CEN"], [2000], "Hydro"], "unitsBuild"] = 0.399  # GW_el
    conv_cap.loc[idx[["HBY"], [2000], "Hydro"], "unitsBuild"] = 0.1422  # GW_el
    conv_cap.loc[idx[["NEL"], [2000], "Hydro"], "unitsBuild"] = 0.0453  # GW_el
    conv_cap.loc[idx[["OTG"], [2000], "Hydro"], "unitsBuild"] = 1.664  # GW_el
    conv_cap.loc[idx[["WTO"], [2000], "Hydro"], "unitsBuild"] = 1.0873  # GW_el
    # do not build anything new
    conv_cap.loc[idx[:, :, :], "noExpansion"] = 1  # boolean
    m.parameter.add(conv_cap, "converter_capacityparam")


    conv_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [hydro_techs, hydro_vintage, hydro_activities, ["Water_in", "Water_out", "Elec"]]
        )
    )
    #FIXME: DIRTY hack TO MATCH 2024 TO REAL DATA MBIE
    conv_coef.loc[idx[:, :, "Power_gen", "Elec"], "coefficient"] = 1 # GW_el
    conv_coef.loc[idx[:, :, "Power_gen", "Water_in"], "coefficient"] = -0.95  # GW_el
    conv_coef.loc[idx[:, :, "Power_gen", "Water_out"], "coefficient"] = 0.95 # GW_el
    #spill is not limited by the capacity of the turbine
    conv_coef.loc[idx[:, :, "Spill", "Water_in"], "coefficient"] = -100 # GW_el
    conv_coef.loc[idx[:, :, "Spill", "Water_out"], "coefficient"] = 100 # GW_el
    # we dont need this bc its just zero: conv_coef.loc[idx[:, :, "Spill", "Elec"], "coefficient"] = 0 # GW_el
    m.parameter.add(conv_coef, "converter_coefficient")

    conv_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], hydro_techs, hydro_vintage]
        )
    )
    conv_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] = 2560 # million EUR / unit
    conv_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1  # binary yes/no
    conv_acc.loc[idx["Invest", :, :, :], "amorTime"] = 20  # years
    conv_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06  # percent/100
    conv_acc.loc[idx["OMFix", :, :, :], "perUnitTotal"] = (
        conv_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] *0.0300
    )  # Mio EUR per unit
    m.parameter.add(conv_acc, "accounting_converterunits")

 
    stor_techs = ["Hydro_reservoir"]

    stor_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([stor_techs, hydro_vintage])
    )
    stor_tech.loc[idx[:, :], "lifeTime"] = 100
    stor_tech.loc[idx[:, :], "levelUpperLimit"] = 1

    m.parameter.add(stor_tech, "storage_techparam")
    stor_tech

	#test nodes
    stor_size = pd.DataFrame(
        index=pd.MultiIndex.from_product([stor_techs,hydro_vintage, ["Water_in"]])
    )
	# question about the units
    stor_size.loc[idx["Hydro_reservoir", :, "Water_in"], "size"] = 1  # GWh_ch / unit  
    stor_size.loc[idx[:, :, "Water_in"], "selfdischarge"] = 0 
    m.parameter.add(stor_size, "storage_sizeparam")

    
    stor_res = pd.DataFrame(
        index=pd.MultiIndex.from_product([hydro_nodes, hydro_years, stor_techs])
    )


    
    # storage reservoir capacity 
    stor_res.loc[idx[:, :, :], "unitsUpperLimit"] = 3000  # units

    stor_res.loc[idx[["CAN"], [2000], :], "unitsBuild"] = 2517.2429 # GWh_el
    stor_res.loc[idx[["HBY"], [2000], :], "unitsBuild"] = 154.2635  # GWh_el
    stor_res.loc[idx[["OTG"], [2000], :], "unitsBuild"] = 729.5595 # GWh_el
    stor_res.loc[idx[["WTO"], [2000], :], "unitsBuild"] = 587.1371 # GWh_el

    stor_res.loc[idx[:, :, :], "noExpansion"] = 1
    m.parameter.add(stor_res, "storage_reservoirparam")

    stor_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], stor_techs, hydro_vintage]
        )
    )
    stor_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] = 1650  # million EUR / unit
    stor_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1  # binary yes/no
    stor_acc.loc[idx["Invest", :, :, :], "amorTime"] = 20  # years
    stor_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06  # percent/100
    stor_acc.loc[idx["OMFix", :, :, :], "perUnitTotal"] = (
        stor_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] * 0.0300
    )  # Mio EUR per unit
    m.parameter.add(stor_acc, "accounting_storageunits")


# conventional

def add_thermal(m):
    the_inst_csv = pd.read_csv(Path(path_brownfield).joinpath("power-plant-nz-database.csv"))
    
    the_vintage = [1950, 2030]
    the_techs = ["BIO", "COAL", "DIE"]     
    the_nodes = [n for n in m.set.nodesdata]    
    the_activities = ["Powergen"]

    the_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([the_techs, the_vintage])
    )
    the_tech.loc[idx[:, :], "lifeTime"] = 30                               # lifeTIME same as gas turbines
    the_tech.loc[idx[:, :], "activityUpperLimit"] = 1                      # UpperLimit same as gas turbines 
    m.parameter.add(the_tech, "converter_techparam") 

    # fix: m.alltheyears or smth for the list - not years calculate
    the_cap = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [the_nodes, range(1992, 2050), the_techs]
        )
    )

    # fixed capactities     
    df = the_inst_csv
    
    filtered_df = df[(df['Type'] == 'Thermal') & (df['Primary_fuel'].isin(['Biogas', 'Biomass', 'Coal', 'Diesel', 'Waste heat', 'Wood', 'Wood waste']))]
    grouped_df = filtered_df.groupby(['Node', 'Year_built', 'Techs'])['Capacity_MW'].sum().reset_index()
    grouped_df['Year_built'] = grouped_df['Year_built'].astype(int) 
    
    the_cap = (grouped_df
               .set_index(["Node", "Year_built", "Techs"])
               .rename(columns={"Capacity_MW": "unitsBuild"})
               .div(1e3))

    the_cap_upper = pd.DataFrame(index=pd.MultiIndex.from_product([the_nodes, yrs_to_calc, the_techs]))
    the_cap_upper.loc[idx[:, [2020], :], "noExpansion"] = 1  # boolean
    the_cap_upper.loc[idx[:, :, :], "unitsUpperLimit"] = 100  # boolean
    the_cap_full = pd.concat([the_cap, the_cap_upper], axis=1)
    m.parameter.add(the_cap_full, "converter_capacityparam")

    
    #MAAM: thermal generation can only convert primary fuel to electricity
    the_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [the_techs, the_vintage, the_activities, ["Elec"]]
        )
    )

    the_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = 1 # GW_el
    m.parameter.add(the_coef, "converter_coefficient")


    the_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], the_techs, the_vintage]
        )
    )
    the_acc.loc[idx["Invest", :, "BIO", :], "perUnitBuild"] = [2600, 2300]  # million EUR / unit
    the_acc.loc[idx["Invest", :, "COAL", :], "perUnitBuild"] = [1600, 1600] # million EUR / unit
    the_acc.loc[idx["Invest", :, "DIE", :], "perUnitBuild"] = [900, 830]    # million EUR / unit
    the_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1                   # binary yes/no
    the_acc.loc[idx["Invest", :, :, :], "amorTime"] = 30                    # years
    the_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06                  # percent/100
    the_acc.loc[idx["OMFix", :, :, :], "perUnitTotal"] = (
        the_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] * 0.033
    )  # Mio EUR per unit
    m.parameter.add(the_acc, "accounting_converterunits")

    # # Emit carbon from combustion
    the_emission = pd.DataFrame(index=pd.MultiIndex.from_product([["CO2_emission", "FuelCost"], ["global"], ["BIO", "COAL","DIE"], the_vintage, the_activities]))
    the_emission.loc[idx["CO2_emission", :, "BIO", :, :], "perActivity"] = np.round(1 / np.array([0.41, 0.43]), 3) * 0.2016 # kt_CO2
    the_emission.loc[idx["CO2_emission", :, "COAL", :, :], "perActivity"] = np.round(1 / np.array([0.58, 0.61]), 3) * 0.2016 # kt_CO2
    the_emission.loc[idx["CO2_emission", :, "DIE", :, :], "perActivity"] = np.round(1 / np.array([0.41, 0.43]), 3) * 0.2016 # kt_CO2
    the_emission.loc[idx["FuelCost", :, "BIO", :, :], "perActivity"] = np.round(1 / np.array([0.41, 0.43]), 3)* 0.2016  # kt_CO2
    the_emission.loc[idx["FuelCost", :, "COAL", :, :], "perActivity"] = np.round(1 / np.array([0.58, 0.61]), 3)* 0.2016 # kt_CO2
    the_emission.loc[idx["FuelCost", :, "DIE", :, :], "perActivity"] = np.round(1 / np.array([0.41, 0.43]), 3)* 0.2016 # kt_CO2
    m.parameter.add(the_emission, "accounting_converteractivity")

def add_gas_turbines(m):    
    gt_inst_csv = pd.read_csv(Path(path_brownfield).joinpath("power-plant-nz-database.csv"))

    gt_vintage = [1950,2030]

    #gt_techs = ["GT", "CCGT", "OCGT", "GT_H2", "CCGT_H2"]                         # OCGT added 
    gt_techs = ["GT", "CCGT", "OCGT"]                       
    gt_nodes = [n for n in m.set.nodesdata if not n.startswith("LNG")]
    gt_activities = ["Powergen"]

    gt_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([gt_techs, gt_vintage])
    )
    gt_tech.loc[idx[:, :], "lifeTime"] = 30
    gt_tech.loc[idx[:, :], "activityUpperLimit"] = 1
    m.parameter.add(gt_tech, "converter_techparam")

    df = gt_inst_csv
    filtered_df = df[df['Primary_fuel'] == 'Natural gas']
    grouped_df = filtered_df.groupby(['Node', 'Year_built', 'Techs'])['Capacity_MW'].sum().reset_index()
    grouped_df['Year_built'] = grouped_df['Year_built'].astype(int) 

    gt_cap = (grouped_df
              .set_index(["Node", "Year_built", "Techs"])
              .rename(columns={"Capacity_MW": "unitsBuild"})
              .div(1e3))

    gt_cap_upper = pd.DataFrame(index=pd.MultiIndex.from_product([gt_nodes, yrs_to_calc, gt_techs]))
    gt_cap_upper.loc[idx[:, :, gt_techs], "unitsUpperLimit"] = 100  # GW_el
    gt_cap_upper.loc[idx[:, [2020], :], "noExpansion"] = 1  # boolean
    gt_cap_full = pd.concat([gt_cap, gt_cap_upper], axis=1)
    m.parameter.add(gt_cap_full, "converter_capacityparam")

    # coefficients  
    gt_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [gt_techs, gt_vintage, gt_activities, ["Elec"]]
        )
    )

    gt_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = 1 # GW_el
    m.parameter.add(gt_coef, "converter_coefficient")
    
    
    # gt_coef = pd.DataFrame(
    #     index=pd.MultiIndex.from_product(
    #         [gt_techs, gt_vintage, gt_activities, ["CH4", "H2", "Elec"]]
    #     )
    # )

    # gt_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = 1 # GW_el
    # gt_coef.loc[idx["GT", :, :, "CH4"], "coefficient"] = np.round(-1 / np.array([0.41, 0.43]), 3) # GW_el #-1/efficiency=coef, rounded to 3, line above output of 1 GW
    # gt_coef.loc[idx["CCGT", :, :, "CH4"], "coefficient"] = np.round(-1 / np.array([0.58, 0.61]), 3) # GW_el
    # gt_coef.loc[idx["OCGT", :, :, "CH4"], "coefficient"] = np.round(-1 / np.array([0.4695, 0.4695]), 3) # GW_el - data: REMix Tutorial 202, number for 2020 used for 2000 and 2030
    # gt_coef.loc[idx["GT_H2", :, :, "H2"], "coefficient"] = np.round(-1 / np.array([.41, 0.43]), 3) # GW_elel
    # gt_coef.loc[idx["CCGT_H2", :, :, "H2"], "coefficient"] = np.round(-1 / np.array([0.58, 0.61]), 3) # GW_el
    # m.parameter.add(gt_coef, "converter_coefficient")

    gt_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], gt_techs, gt_vintage]
        )
    )
    gt_acc.loc[idx["Invest", :, "GT", :], "perUnitBuild"] = [900, 830]  # million EUR / unit
    gt_acc.loc[idx["Invest", :, "CCGT", :], "perUnitBuild"] = [775, 775]  # million EUR / unit
    gt_acc.loc[idx["Invest", :, "OCGT", :], "perUnitBuild"] = [475, 475]  # million EUR / unit  - data: LUT Breyer "Financial_Technical assumptions-newversion.docx"
    #gt_acc.loc[idx["Invest", :, "GT_H2", :], "perUnitBuild"] = [900, 830]  # million EUR / unit
    #gt_acc.loc[idx["Invest", :, "CCGT_H2", :], "perUnitBuild"] = [600, 560]  # million EUR / unit
    gt_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1  # binary yes/no
    gt_acc.loc[idx["Invest", :, :, :], "amorTime"] = 30  # years
    gt_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06  # percent/100
    gt_acc.loc[idx["OMFix", :, :, :], "perUnitTotal"] = (
        gt_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] * 0.0193
    )  # Mio EUR per unit
    m.parameter.add(gt_acc, "accounting_converterunits")

    # # Emit carbon from combustion
    gt_emission = pd.DataFrame(index=pd.MultiIndex.from_product([["CO2_emission", "FuelCost"], ["global"], ["GT", "CCGT","OCGT"], gt_vintage, gt_activities]))
    gt_emission.loc[idx["CO2_emission", :, "GT", :, :], "perActivity"] = np.round(1 / np.array([0.41, 0.43]), 3) * 0.2016 # kt_CO2
    gt_emission.loc[idx["CO2_emission", :, "CCGT", :, :], "perActivity"] = np.round(1 / np.array([0.58, 0.61]), 3) * 0.2016 # kt_CO2
    gt_emission.loc[idx["CO2_emission", :, "OCGT", :, :], "perActivity"] = np.round(1 / np.array([0.4695, 0.4695]), 3) * 0.2016 # kt_CO2 - data: REMix Tutorial 202, number for 2020 used for 2000 and 2030
    gt_emission.loc[idx["FuelCost", :, "GT", :, :], "perActivity"] = np.round(1 / np.array([0.41, 0.43]), 3)* 0.2016 # kt_CO2
    gt_emission.loc[idx["FuelCost", :, "CCGT", :, :], "perActivity"] = np.round(1 / np.array([0.58, 0.61]), 3)* 0.2016 # kt_CO2
    gt_emission.loc[idx["FuelCost", :, "OCGT", :, :], "perActivity"] = np.round(1 / np.array([0.4695, 0.4695]), 3)* 0.2016 # kt_CO2 - data: REMix Tutorial 202, number for 2020 used for 2000 and 2030
    # Emit carbon from combustion
    # gt_emission = pd.DataFrame(index=pd.MultiIndex.from_product([["CO2_emission", "FuelCost"], ["global"], ["GT", "CCGT","OCGT"], gt_vintage, gt_activities]))
    # gt_emission.loc[idx["CO2_emission", :, "GT", :, :], "perActivity"] = np.round(1 / np.array([0.41]), 3) * 0.2016 # kt_CO2
    # gt_emission.loc[idx["CO2_emission", :, "CCGT", :, :], "perActivity"] = np.round(1 / np.array([0.58]), 3) * 0.2016 # kt_CO2
    # gt_emission.loc[idx["CO2_emission", :, "OCGT", :, :], "perActivity"] = np.round(1 / np.array([0.4695]), 3) * 0.2016 # kt_CO2 - data: REMix Tutorial 202, number for 2020 used for 2000 and 2030
    # gt_emission.loc[idx["FuelCost", :, "GT", :, :], "perActivity"] = np.round(1 / np.array([0.41]), 3)* 0.2016 # kt_CO2
    # gt_emission.loc[idx["FuelCost", :, "CCGT", :, :], "perActivity"] = np.round(1 / np.array([0.58]), 3)* 0.2016 # kt_CO2
    # gt_emission.loc[idx["FuelCost", :, "OCGT", :, :], "perActivity"] = np.round(1 / np.array([0.4695]), 3)* 0.2016 # kt_CO2 - data: REMix Tutorial 202, number for 2020 used for 2000 and 2030
    
    m.parameter.add(gt_emission, "accounting_converteractivity")

# hydrogen

def add_electrolyser(m):
    eltr_vintage = [2020, 2030, 2040, 2050]
    eltr_nodes = [n for n in m.set.nodesdata]
    # technology
    eltr_tech = pd.DataFrame(index=pd.MultiIndex.from_product([["Electrolyser"], eltr_vintage]))
    eltr_tech.loc[idx[:, :], ["lifeTime"]] = [25, 30, 32, 35]  # years
    eltr_tech["activityUpperLimit"] = 1  # availability of technology
    m.parameter.add(eltr_tech, "converter_techparam")

    # capacities
    eltr_caps = pd.DataFrame(index=pd.MultiIndex.from_product([eltr_nodes, yrs_to_calc, ["Electrolyser"]]))
    eltr_caps["unitsUpperLimit"] = 100  # GW_el
    m.parameter.add(eltr_caps, "converter_capacityparam")

    # coefficients
    eltr_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [
                ["Electrolyser"],
                eltr_vintage,
                ["Electrolysis"],
                ["Elec", "H2"],
            ]
        )
    )
    eltr_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = -1
    eltr_coef.loc[idx[:, :, :, "H2"], "coefficient"] = [0.665, 0.79, 0.79, 0.85] # DEA2022 AEC 1MW comm&indust for 2020, Will for the others
    m.parameter.add(eltr_coef, "converter_coefficient")

    # accounting
    electrolyser_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product([["Invest", "OMFix"], ["global"], ["Electrolyser"], eltr_vintage])
        ).sort_index()

    electrolyser_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] = [750, 570, 450, 350]   #  DEA2022 AEC 1MW comm&indust - Mio EUR per unit
    electrolyser_acc.loc[idx["Invest", :, :, :], "amorTime"] = [25, 30, 32, 35]  # years
    electrolyser_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1  # binary yes/no
    electrolyser_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06  # percent/100
    electrolyser_acc.loc[idx["OMFix", :, :, :], "perUnitTotal"] = (
        electrolyser_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] * 0.014
        )  # Mio EUR per unit

    m.parameter.add(electrolyser_acc, "accounting_converterunits")


def add_dac(m):
    dac_vintage = [2020, 2030, 2040, 2050]
    dac_techs = ["DAC"]
    dac_nodes = [n for n in m.set.nodesdata if not n.startswith("LNG")]

    # technology
    dac_tech = pd.DataFrame(index=pd.MultiIndex.from_product([dac_techs, dac_vintage]))
    dac_tech.loc[idx[:, [2020, 2030, 2040, 2050]], ["lifeTime"]] = [20]  # years
    dac_tech["activityUpperLimit"] = 1  # availability of technology
    m.parameter.add(dac_tech, "converter_techparam")

    # capacities
    dac_caps = pd.DataFrame(index=pd.MultiIndex.from_product([dac_nodes, [2020], dac_techs]))
    dac_caps["unitsUpperLimit"] = 100  # GW_el
    m.parameter.add(dac_caps, "converter_capacityparam")

    # coefficients
    dac_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [
                dac_techs,
                dac_vintage,
                ["Capture"],
                ["Elec", "CO2_feed"],
            ]
        )
    )
    dac_coef.loc[idx[:, :, :, "CO2_feed"], "coefficient"] = 1
    dac_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = [-1.535, -1.458, -1.385, -1.316] # GWh/el per ktCO2
    m.parameter.add(dac_coef, "converter_coefficient")

    # accounting
    dac_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product([["Invest", "OMFix"], ["global"], dac_techs, dac_vintage])
    ).sort_index()

    dac_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] = [i * 8.76 for i in [815]]  # EUR/tCO2*a -> MEUR/ktCO2*h
    dac_acc.loc[idx["Invest", :, :, :], "amorTime"] = [20]  # years
    dac_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1  # binary yes/no
    dac_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06  # percent/100
    dac_acc.loc[idx["OMFix", :, :, :], "perUnitTotal"] = (
        dac_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] * 0.04
    )  # Mio EUR per unit

    m.parameter.add(dac_acc, "accounting_converterunits")

    # Remove Carbon from CO2_emission indicator
    dac_activity = pd.DataFrame(index=pd.MultiIndex.from_product([["CO2_emission"], ["global"], ["DAC"], dac_vintage, ["Capture"]]))
    dac_activity["perActivity"] = -1
    m.parameter.add(dac_activity, "accounting_converteractivity")

def add_methanizer(m):
    methanizer_vintage = [2020, 2030, 2040, 2050]
    methanizer_techs = ["Methanizer"]
    methanizer_nodes = [n for n in m.set.nodesdata if not n.startswith("LNG")]

    # technology
    methanizer_tech = pd.DataFrame(index=pd.MultiIndex.from_product([methanizer_techs, methanizer_vintage]))
    methanizer_tech.loc[idx[:, :], ["lifeTime"]] = 30  # years
    methanizer_tech["activityUpperLimit"] = 1  # availability of technology
    m.parameter.add(methanizer_tech, "converter_techparam")

    # capacities
    methanizer_caps = pd.DataFrame(index=pd.MultiIndex.from_product([methanizer_nodes, [2020, 2025, 2030, 2035, 2040, 2045, 2050], methanizer_techs]))
    methanizer_caps["unitsUpperLimit"] = 100  # GW_el
    m.parameter.add(methanizer_caps, "converter_capacityparam")

    # coefficients
    methanizer_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [
                methanizer_techs,
                methanizer_vintage,
                methanizer_techs,
                ["Elec", "H2", "CH4", "CO2_feed"],
            ]
        )
    )
    methanizer_coef.loc[idx[:, :, :, "CH4"], "coefficient"] = 1
    methanizer_coef.loc[idx[:, :, :, "H2"], "coefficient"] = -1.284
    methanizer_coef.loc[idx[:, :, :, "CO2_feed"], "coefficient"] = -0.2016
    methanizer_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = -0.006

    m.parameter.add(methanizer_coef, "converter_coefficient")

    # accounting
    methanizer_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product([["Invest", "OMFix"], ["global"], methanizer_techs, methanizer_vintage])
    ).sort_index()

    methanizer_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] = [558, 309, 251, 211]  # Mio EUR per unit
    methanizer_acc.loc[idx["Invest", :, :, :], "amorTime"] = [20, 25, 30, 30]  # years
    methanizer_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1  # binary yes/no
    methanizer_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06  # percent/100
    methanizer_acc.loc[idx["OMFix", :, :, :], "perUnitTotal"] = (
        methanizer_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] * 0.04
    )  # Mio EUR per unit

    m.parameter.add(methanizer_acc, "accounting_converterunits")

def add_methanol_syn(m):
    methanol_syn_vintage = [2020]
    methanol_syn_techs = ["MethanolSyn"]
    methanol_nodes = [n for n in m.set.nodesdata if not n.startswith("LNG")]

    # technology
    methanol_syn_tech = pd.DataFrame(index=pd.MultiIndex.from_product([methanol_syn_techs, methanol_syn_vintage]))
    methanol_syn_tech.loc[idx[:, :], ["lifeTime"]] = 30  # years
    methanol_syn_tech["activityUpperLimit"] = 1  # availability of technology
    m.parameter.add(methanol_syn_tech, "converter_techparam")

    # capacities
    methanol_syn_caps = pd.DataFrame(index=pd.MultiIndex.from_product([methanol_nodes, [2020, 2025, 2030, 2035, 2040, 2045, 2050], methanol_syn_techs]))
    methanol_syn_caps["unitsUpperLimit"] = 100  # GW_el
    m.parameter.add(methanol_syn_caps, "converter_capacityparam")

    # coefficients
    methanol_syn_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [
                methanol_syn_techs,
                methanol_syn_vintage,
                ["Synthesis"],
                ["CH3OH", "H2", "CO2_feed", "Elec"],
            ]
        )
    )
    methanol_syn_coef.loc[idx[:, :, :, "CH3OH"], "coefficient"] = 1
    methanol_syn_coef.loc[idx[:, :, :, "H2"], "coefficient"] = -1.25
    methanol_syn_coef.loc[idx[:, :, :, "CO2_feed"], "coefficient"] = -0.219
    methanol_syn_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = -0.1

    m.parameter.add(methanol_syn_coef, "converter_coefficient")

    # accounting
    methanol_syn_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product([["Invest", "OMFix"], ["global"], methanol_syn_techs, methanol_syn_vintage])
    ).sort_index()

    methanol_syn_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] = 835  # Mio EUR per unit
    methanol_syn_acc.loc[idx["Invest", :, :, :], "amorTime"] = 30  # years
    methanol_syn_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1  # binary yes/no
    methanol_syn_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06  # percent/100
    methanol_syn_acc.loc[idx["OMFix", :, :, :], "perUnitTotal"] = (
        methanol_syn_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] * 0.04
    )  # Mio EUR per unit

    m.parameter.add(methanol_syn_acc, "accounting_converterunits")

def add_ftropsch_syn(m):
    # TODO: Produce syncrude for refinieries 
    ftropsch_vintage = [2020, 2040]
    ftropsch_techs = ["FTropschSyn"]
    ftropsch_nodes = [n for n in m.set.nodesdata if not n.startswith("LNG")]

    # technology
    ftropsch_syn_tech = pd.DataFrame(index=pd.MultiIndex.from_product([ftropsch_techs, ftropsch_vintage]))
    ftropsch_syn_tech.loc[idx[:, :], ["lifeTime"]] = [30, 30]  # years
    ftropsch_syn_tech["activityUpperLimit"] = 1  # availability of technology
    m.parameter.add(ftropsch_syn_tech, "converter_techparam")

    # capacities
    ftropsch_syn_caps = pd.DataFrame(index=pd.MultiIndex.from_product([ftropsch_nodes, [2020, 2025, 2030, 2035, 2040, 2045, 2050], ftropsch_techs]))
    ftropsch_syn_caps["unitsUpperLimit"] = 100  # GW_el
    m.parameter.add(ftropsch_syn_caps, "converter_capacityparam")

    # coefficients - data from Andi SynLink? 
    ftropsch_syn_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [
                ftropsch_techs,
                ftropsch_vintage,
                ["Synthesis"],
                ["REfuel", "H2", "CO2_feed", "Elec"], # "H2O"
            ]
        )
    )
    # Esitmates Andi, 65% synthesis efficiency
    ftropsch_syn_coef.loc[idx[:, :, :, "REfuel"], "coefficient"] = 1
    ftropsch_syn_coef.loc[idx[:, :, :, "H2"], "coefficient"] = -1.52
    ftropsch_syn_coef.loc[idx[:, :, :, "CO2_feed"], "coefficient"] = -0.35
    ftropsch_syn_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = -0.11
    # ftropsch_syn_coef.loc[idx[:, :, :, "H2O"], "coefficient"] = 1

    m.parameter.add(ftropsch_syn_coef, "converter_coefficient")

    # accounting
    ftropsch_syn_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product([["Invest", "OMFix"], ["global"], ftropsch_techs, ftropsch_vintage])
    ).sort_index()

    ftropsch_syn_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] = [1017, 915]  # Mio EUR per unit
    ftropsch_syn_acc.loc[idx["Invest", :, :, :], "amorTime"] = [30, 30]  # years
    ftropsch_syn_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1  # binary yes/no
    ftropsch_syn_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06  # percent/100
    ftropsch_syn_acc.loc[idx["OMFix", :, :, :], "perUnitTotal"] = (
        ftropsch_syn_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] * 0.04
    )  # Mio EUR per unit

    m.parameter.add(ftropsch_syn_acc, "accounting_converterunits")

# storage
    
def add_lithium_batteries(m):
    battery_vintage = [2020, 2030, 2040, 2050] #FIXME: 2035 was 2030, change when not running will's values
    battery_techs = ["Battery"]
    battery_nodes = [n for n in m.set.nodesdata if not n.startswith("LNG")]

    conv_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([battery_techs, battery_vintage])
    )
    conv_tech.loc[idx[:, :], "lifeTime"] = 20
    conv_tech.loc[idx[:, :], "activityUpperLimit"] = 1
    m.parameter.add(conv_tech, "converter_techparam")

    conv_cap = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [battery_nodes, yrs_sel, battery_techs]
        )
    )
    conv_cap.loc[idx[:, :, battery_techs], "unitsUpperLimit"] = 50  # GW_el
    m.parameter.add(conv_cap, "converter_capacityparam")


    conv_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [battery_techs,battery_vintage, ["Charge", "Discharge"], ["Elec", "Elec_battery"]]
        )
    )
    conv_coef.loc[idx[:, :, "Charge", "Elec"], "coefficient"] = -1 # GW_el
    conv_coef.loc[idx[:, :, "Charge", "Elec_battery"], "coefficient"] = 0.975 # GW_el
    conv_coef.loc[idx[:, :, "Discharge", "Elec"], "coefficient"] = 1 # GW_el
    conv_coef.loc[idx[:, :, "Discharge", "Elec_battery"], "coefficient"] = -1.025 # GW_el
    m.parameter.add(conv_coef, "converter_coefficient")

    conv_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], battery_techs, battery_vintage]
        )
    )
    conv_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] = [117, 55, 37, 30]  # million EUR / unit
    conv_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1  # binary yes/no
    conv_acc.loc[idx["Invest", :, :, :], "amorTime"] = 20  # years
    conv_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06  # percent/100
    conv_acc.loc[idx["OMFix", :, :, :], "perUnitTotal"] = (
        conv_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] * 0.014
    )  # Mio EUR per unit
    m.parameter.add(conv_acc, "accounting_converterunits")


    stor_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([battery_techs, battery_vintage])
    )
    stor_tech.loc[idx[:, :], "lifeTime"] = 20
    stor_tech.loc[idx[:, :], "levelUpperLimit"] = 1

    m.parameter.add(stor_tech, "storage_techparam")
    stor_tech


    stor_size = pd.DataFrame(
        index=pd.MultiIndex.from_product([battery_techs, battery_vintage, ["Elec_battery"]])
    )
    stor_size.loc[idx["Battery", :, "Elec_battery"], "size"] = 4  # GWh_ch / unit
    m.parameter.add(stor_size, "storage_sizeparam")


    stor_res = pd.DataFrame(
        index=pd.MultiIndex.from_product([battery_nodes,  yrs_sel, battery_techs])
    )
    stor_res.loc[idx[:, :, :], "unitsUpperLimit"] = 30  # units
    stor_res.loc[idx[:, [2020], :], "noExpansion"] = 1  # boolean
    m.parameter.add(stor_res, "storage_reservoirparam")

    stor_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], battery_techs, battery_vintage]
        )
    )
    stor_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] = [i * 4 for i in [234, 110, 76, 61]]  # million EUR / unit
    stor_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1  # binary yes/no
    stor_acc.loc[idx["Invest", :, :, :], "amorTime"] = 20  # years
    stor_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06  # percent/100
    stor_acc.loc[idx["OMFix", :, :, :], "perUnitTotal"] = (
        stor_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] * 0.014
    )  # Mio EUR per unit
    m.parameter.add(stor_acc, "accounting_storageunits")

def add_h2_storage(m):

    # converter is the compressor 
    # storage is the gas tank
    h2_stor_vintage = [2020]#, 2030, 2040, 2050]
    yrs_h2=[2020]#, 2025, 2030, 2035, 2040, 2045, 2050]
    h2_stor_techs = ["H2_storage"]
    h2_stor_nodes = [n for n in m.set.nodesdata]

    conv_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([h2_stor_techs, yrs_h2])
    )
    conv_tech.loc[idx[:, :], "lifeTime"] = 40
    conv_tech.loc[idx[:, :], "activityUpperLimit"] = 1
    m.parameter.add(conv_tech, "converter_techparam")

    conv_cap = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [h2_stor_nodes, yrs_h2, h2_stor_techs]
        )
    )
    conv_cap.loc[idx[:, :, h2_stor_techs], "unitsUpperLimit"] = 50 # GW_el
    m.parameter.add(conv_cap, "converter_capacityparam")


    conv_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [h2_stor_techs,  yrs_h2, ["Charge", "Discharge"], ["H2", "H2_stored","Elec"]]
        )
    )
    # electrolyzer_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = -1
    conv_coef.loc[idx[:, :, "Charge", "Elec"], "coefficient"] = [ -0.043346085]#,	-0.043346085, -0.035470537,	-0.035470537, -0.035470537, -0.031529343, -0.031529343]
    # 2020, 2025, 030, 2035, 2040, 2045, 2040...
    conv_coef.loc[idx[:, :, "Charge", "H2"], "coefficient"] = -0.15  # GW_h2
    conv_coef.loc[idx[:, :, "Charge", "H2_stored"], "coefficient"] = 0.15  * 0.89 # GW_h2
    conv_coef.loc[idx[:, :, "Discharge", "H2"], "coefficient"] = 0.15  # GW_h2
    conv_coef.loc[idx[:, :, "Discharge", "H2_stored"], "coefficient"] = -0.15 * 1.11 # GW_h2
    m.parameter.add(conv_coef, "converter_coefficient")

    conv_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], h2_stor_techs, yrs_h2]
        )
    )
    #
    conv_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] = 5.14  # million EUR / unit
    conv_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1  # binary yes/no
    conv_acc.loc[idx["Invest", :, :, :], "amorTime"] = 40  # years
    conv_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06  # percent/100
    conv_acc.loc[idx["OMFix", :, :, :], "perUnitTotal"] = (
        conv_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] * 0.04
    )  # Mio EUR per unit
    m.parameter.add(conv_acc, "accounting_converterunits")


    stor_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([h2_stor_techs, yrs_h2])
    )
    stor_tech.loc[idx[:, :], "lifeTime"] = 40
    stor_tech.loc[idx[:, :], "levelUpperLimit"] = 1

    m.parameter.add(stor_tech, "storage_techparam")
    stor_tech


    stor_size = pd.DataFrame(
        index=pd.MultiIndex.from_product([h2_stor_techs, yrs_h2, ["H2_stored"]])
    )
    stor_size.loc[idx["H2_storage", :, "H2_stored"], "size"] = 6.9993  # GWh / unit
    m.parameter.add(stor_size, "storage_sizeparam")


    stor_res = pd.DataFrame(
        index=pd.MultiIndex.from_product([h2_stor_nodes, yrs_h2, h2_stor_techs])
    )
    stor_res.loc[idx[:, :, :], "unitsUpperLimit"] = 50  # units 
    stor_res.loc[idx[:, [2020], :], "noExpansion"] = 1  # boolean
    m.parameter.add(stor_res, "storage_reservoirparam")

    stor_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], h2_stor_techs, yrs_h2]
        )
    )
   
     #CAPEX (628.14 EUR / kgH2 * 210000 kgH2/unit ) /1M = 131.91 M EUR / unit
    stor_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] = 131.91  # million EUR / unit
    stor_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1  # binary yes/no
    stor_acc.loc[idx["Invest", :, :, :], "amorTime"] = 40  # years 
    stor_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06  # percent/100
    stor_acc.loc[idx["OMFix", :, :, :], "perUnitTotal"] = (
        stor_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] * 0.03
    )  # Mio EUR per unit
    m.parameter.add(stor_acc, "accounting_storageunits")

# others
    
def add_network(m):
    #First we need to set up the link connections in the data by defining the starting and ending node of each link
    link_names = ["NIS__AKL",
                  "AKL__WTO",
                  "WTO__BOP",
                  "WTO__CEN",
                  "CEN__HBY",
                  "TRN__CEN",
                  "CEN__WEL",
                  "WEL__CAN",
                  "NEL__CAN",
                  "CAN__OTG",
                  "AKL__TRN",
                  "WTO__HBY"]
    nodes_data = set(m.set.nodesdata)
    
    link_connections = pd.DataFrame(index=pd.MultiIndex.from_product([link_names, nodes_data]))
    link_connections.loc[idx["NIS__AKL", "NIS"], ["start"]] = 1
    link_connections.loc[idx["NIS__AKL", "AKL"], ["end"]] = 1
    link_connections.loc[idx["AKL__WTO","AKL"], ["start"]] = 1
    link_connections.loc[idx["AKL__WTO", "WTO"], ["end"]] = 1
    link_connections.loc[idx["WTO__BOP", "WTO"], ["start"]] = 1
    link_connections.loc[idx["WTO__BOP", "BOP"], ["end"]] = 1
    link_connections.loc[idx["WTO__CEN", "WTO"], ["start"]] = 1
    link_connections.loc[idx["WTO__CEN", "CEN"], ["end"]] = 1
    link_connections.loc[idx["CEN__HBY", "CEN"], ["start"]] = 1
    link_connections.loc[idx["CEN__HBY", "HBY"], ["end"]] = 1
    link_connections.loc[idx["TRN__CEN", "TRN"], ["start"]] = 1
    link_connections.loc[idx["TRN__CEN", "CEN"], ["end"]] = 1
    link_connections.loc[idx["CEN__WEL", "CEN"], ["start"]] = 1
    link_connections.loc[idx["CEN__WEL", "WEL"], ["end"]] = 1
    link_connections.loc[idx["WEL__CAN", "WEL"], ["start"]] = 1
    link_connections.loc[idx["WEL__CAN", "CAN"], ["end"]] = 1
    link_connections.loc[idx["NEL__CAN", "NEL"], ["start"]] = 1
    link_connections.loc[idx["NEL__CAN", "CAN"], ["end"]] = 1
    link_connections.loc[idx["CAN__OTG", "CAN"], ["start"]] = 1
    link_connections.loc[idx["CAN__OTG", "OTG"], ["end"]] = 1
    link_connections.loc[idx["AKL__TRN", "AKL"], ["start"]] = 1
    link_connections.loc[idx["AKL__TRN", "TRN"], ["end"]] = 1
    link_connections.loc[idx["WTO__HBY", "WTO"], ["start"]] = 1
    link_connections.loc[idx["WTO__HBY", "HBY"], ["end"]] = 1
    link_connections = link_connections.dropna(how="all").fillna(0)

    m.parameter.add(link_connections, "transfer_linkstartend")
    link_connections

    #Next we define the lengths for each corridor. We can use different distance types.
    # "transport_distanceParam"
    link_types = ["land", "sea"]
 
    link_lengths = pd.DataFrame(index=pd.MultiIndex.from_product([link_names, link_types]))
    link_lengths.loc[idx["NIS__AKL", "land"], ["length"]] = 149.8
    link_lengths.loc[idx["AKL__WTO", "land"], ["length"]] = 136.0
    link_lengths.loc[idx["WTO__BOP", "land"], ["length"]] = 76.2
    link_lengths.loc[idx["WTO__CEN", "land"], ["length"]] = 154.5
    link_lengths.loc[idx["CEN__HBY", "land"], ["length"]] = 96.5
    link_lengths.loc[idx["TRN__CEN", "land"], ["length"]] = 115.1
    link_lengths.loc[idx["CEN__WEL", "land"], ["length"]] = 83.1
    link_lengths.loc[idx["WEL__CAN", "land"], ["length"]] = 316.4
    link_lengths.loc[idx["NEL__CAN", "land"], ["length"]] = 203.9
    link_lengths.loc[idx["CAN__OTG", "land"], ["length"]] = 179.5
    link_lengths.loc[idx["AKL__TRN", "land"], ["length"]] = 200.1
    link_lengths.loc[idx["WTO__HBY", "land"], ["length"]] = 114.6
    # link_lengths = link_lengths.dropna()
    m.parameter.add(link_lengths, "transfer_lengthparam")

    transfer_lengths = pd.DataFrame(
    index=pd.MultiIndex.from_product([link_names, link_types]))

    #With the line corridors now defined, we can start adding lines to be optimized to the model.
    #Note: we have to change these
    transport_techs = ["HV"] #other options: "HVDC", "Pipeline_CH4", "Pipeline_H2", "Pipeline_H2_retrofit"
    #in previous version of this function, the capacities come from a csv file
    # question: should I have a different capacity for each tech? how do I define that?
    # question: here I have transport_techs, how do I assign different transport techs to different links (to
    # where i want to have the hvDC? also should i just split the line like that and create extra virtual nodes?
    #  or how do i go about that line having only one fragment being DC)
    link_caps = pd.DataFrame(
        index=pd.MultiIndex.from_product([link_names, m.set.yearssel, transport_techs])
    )
    link_caps.loc[idx[:,:,"HV"],["linksUpperLimit"]] = 100 # Allow to build 100 GW for all links as the upper limit
    

    m.parameter.add(link_caps, "transfer_linksparam")
    tech_params = pd.DataFrame(
        index=pd.MultiIndex.from_product([transport_techs, m.set.yearssel])
    )
    tech_params.loc[:, "lifeTime"] = 40
    tech_params.loc[:, "flowUpperLimit"] = 1

    m.parameter.add(tech_params, "transfer_techparam")

    # Define the commodity and rated capacity of the network technology
    # "transfer_coefficient"
    commodity = ["Elec"]

    transfer_coefficient = pd.DataFrame(
        index=pd.MultiIndex.from_product([transport_techs, m.set.yearssel, commodity])
    )
    transfer_coefficient["coefficient"] = 1  # GWh / h per line

    m.parameter.add(transfer_coefficient, "transfer_coefficient")
    transfer_coefficient


    
    # Define the losses for the converter stations
    # "transport_coefPerFlow"
    coef_per_flow = pd.DataFrame(
        index=pd.MultiIndex.from_product([transport_techs, m.set.yearssel, commodity])
    )
    coef_per_flow[
        "coefPerFlow"
    ] = -0.014  # electrical losses of 14 MWh/h for each flow of 1 GWh/h

    m.parameter.add(coef_per_flow, "transfer_coefperflow")

    # Define the losses for the lines per km
    # "transport_coefPerDistance"
    coef_per_dist = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [transport_techs, m.set.yearssel, commodity, link_types]
        )
    )
    coef_per_dist.loc[
        idx[:, :, :, "land"], idx["coefPerLength"]
    ] = (
        -0.00004
    )  # electrical losses of 40 kWh / h for each flow of 1 GWh / h and 1 km line length ~ 24 MWh / h for 600 km distance
    coef_per_dist.loc[idx[:, :, :, "sea"], idx["coefPerLength"]] = -0.00003

    m.parameter.add(coef_per_dist, "transfer_coefperlength")
    coef_per_dist

    # Define indicators for each line built (for HVDC this is an AC/DC converter station 
    # at the beginning and end of the line)
    # "accounting_transferlinks"
    cost_indicators = ["Invest", "OMFix"]
    area = ["global"]

    transfer_indicators = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [cost_indicators, area, transport_techs, m.set.yearssel]
        )
    )
    transfer_indicators.loc[idx["Invest", "global"], "perLinkBuild"] = 180
    transfer_indicators.loc[idx["Invest", "global"], "interest"] = 0.06
    transfer_indicators.loc[idx["Invest", "global"], "amorTime"] = 40
    transfer_indicators.loc[idx["Invest", "global"], "useAnnuity"] = 1
    transfer_indicators.loc[idx["OMFix", "global"], "perLinkTotal"] = 1.8
    transfer_indicators = transfer_indicators.fillna(0)

    m.parameter.add(transfer_indicators, "accounting_transferlinks")
    transfer_indicators

    # Define indicators for each line-km built 
    # (this needs the additional set for distance-type modifiers, such as land and sea)
    # ""accounting_transferperlength""
    indicators_distance = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [cost_indicators, area, transport_techs, m.set.yearssel, link_types]
        )
    )
    indicators_distance.loc[
        idx["Invest", "global", :, :, "land"], "perLengthBuild"
    ] = 0.544
    indicators_distance.loc[idx["Invest", "global", :, :, "land"], "interest"] = 0.06
    indicators_distance.loc[idx["Invest", "global", :, :, "land"], "amorTime"] = 40
    indicators_distance.loc[idx["Invest", "global", :, :, "land"], "useAnnuity"] = 1
    indicators_distance.loc[
        idx["OMFix", "global", :, :, "land"], "perLengthTotal"
    ] = 0.00544

    indicators_distance.loc[
        idx["Invest", "global", :, :, "sea"], "perLengthBuild"
    ] = 0.975
    indicators_distance.loc[idx["Invest", "global", :, :, "sea"], "interest"] = 0.06
    indicators_distance.loc[idx["Invest", "global", :, :, "sea"], "amorTime"] = 40
    indicators_distance.loc[idx["Invest", "global", :, :, "sea"], "useAnnuity"] = 1
    indicators_distance.loc[
        idx["OMFix", "global", :, :, "sea"], "perLengthTotal"
    ] = 0.00975
    indicators_distance = indicators_distance.fillna(0)

    m.parameter.add(indicators_distance, "accounting_transferperlength")
    indicators_distance

def add_accounting(m):
    #  The value global uses all the regions in the system
    #  whereas the value horizon takes into account all years in the set set.yearssel 
    accounting_indicatorBounds = pd.DataFrame(
        index=pd.MultiIndex.from_product([["global"], ["horizon"], ["SystemCost"]])
    )
    accounting_indicatorBounds["obj"] = -1  # minimization of system costs
    accounting_indicatorBounds["discount"] = 0.02  # social discount rate for the indicators
    m.parameter.add(accounting_indicatorBounds, "accounting_indicatorbounds")

    accounting_perIndicator = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [
                ["SystemCost"],
                ["Invest", "OMFix", "FuelCost", "ImportCost", "SlackCost"],
                ["global"],
                ["horizon"],
            ]
        )
    )
    accounting_perIndicator["perIndicator"] = 1
    m.parameter.add(accounting_perIndicator, "accounting_perindicator")

def validate_scope(m):
    m.infer_set_data()
    nodes_data = set(m.set.nodesdata)
    nodes_model = set(m.set.nodesmodel)
    print(f"Not including modes nodes: {', '.join(sorted(nodes_data - nodes_model))}")
s1 = time.perf_counter()
# #%%
if __name__ == "__main__":
    # Create instance
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
    m.write(output_path=data_dir, fileformat="csv")
    e2 = time.perf_counter()
    d2=time.strftime("%Hh %Mm %Ss", time.gmtime(e2-s2))
    print(f"Writing dataset for {case_name} took {d2}.")
# %%
    # m_run = Instance.from_path(data_dir)

    #COMMENT: YOU NEED TO  COMMENT THE LINE ABOVE SCEN IS DEFINED HERE AND THEN THE TWO LINES BELOW ARE ALSO CHANGED
    #scen = "Wind+"

    m.run(
        resultdir=results_dir,
        resultfile=f"{case_name}",
        #resultfile=f"{case_name}_{scen.replace('/', '_')}",
        #scendir=scen,
        solver="gurobi",
        threads=8,
        lo=4,
        timeres=1,
        names=1,
        roundts=1,
        iis=1,
        # profile=1,
        # gdx="default",
        pathopt="myopic",
        barorder=1
    )
    print(os. getcwd())

e1 = time.perf_counter()
d1=time.strftime("%Hh %Mm %Ss", time.gmtime(e1-s1))
print(f"------------- Running {case_name} took {d1}.")
