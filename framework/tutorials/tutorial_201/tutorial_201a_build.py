# %% [markdown]
# (tutorial_201_label)=
#
# # Tutorial 201 - Electric vehicles
#
# <div style="text-align: center;">
#
# ![Model overview for tutorial 201](../../img/REMix_tutorial201.svg "Model overview for tutorial 201")
#
# Model overview per region of tutorial 201
#
# </div>
#
# ## Part a: setting up the model
#
# In this tutorial we have a closer look at **electric vehicles** with
# controlled and uncontrolled charging.
#
# We will use the model from tutorial_103 as a base model by reading its files
# into a Instance object `m` and adding electric vehicles (EVs) to it in every
# model node.

# %%
# importing dependencies
from remix.framework import Instance
import pandas as pd
import numpy as np
import pathlib as pt

# reading tutorial 103 into REMixInstance object `m`
_path_tut3_data = pt.Path("../tutorial_103/data")

if not _path_tut3_data.exists():
    raise IOError("You need to run tutorial 103a first!")

m = Instance.from_path(_path_tut3_data)

m.datadir = "./data"

# define often-used shortcut
idx = pd.IndexSlice
# %% [markdown]
# In this tutorial, we will both model controlled and uncontrolled charging of the vehicles.
# Controlled and uncontrolled charging make use of the storage feature in the framework. This means, the vehicle
# fleet is modeled as a stationary battery storage system (including a storage reservoir and converters for charging
# and discharging) and a demand profile (a fixed demand that always needs to be satisfied, which represent the actual
# electricity needed for driving purposes).

# %% [markdown]
# EV Storage Converter (Charging Infrastucture - Energy System) for Controlled Charging (CC)

# Power:
# uncontrolled charging  = values for an average vehicle - unit: kW/vehicle (-> need vehicle number per node)
# drivePower = values for an average vehicle - unit: kW/vehicle (-> need vehicle number per node)
# maxChargeAvail = average maximum hourly recharge capacity of the EV fleet - unit: kW/vehicle
# (given in fleet percentage -> need vehicle number per node)

# Energy:
# State profiles (SoC mon and SoC max) = profiles containing the minimum energy level of the vehicle batteries.
# Profiles can be either in energy units [kWh] or relative to the battery capacity.
# Activity profile in REMix are percentages [0:1] of maximum converter capacity
# Coefficient profile (dynamic coefficients)
# Approach: converter vehicle base into capacity base
# (potentially include this transformation into VencoPy), remove vehicle side from REMix, rather a VencoPy thing

# Profiles are normalized, need additional inputs related to EV annual demand to scale and
# share of controlled charging (CC) & uncontrolled charging (UC)
# %% [markdown]
# First we define the electricity demand of the overall car park of the corresponding vehicle technology
# and region [TWh/a] and the vehicle number per model node (which ideally differs per model year).
# %%
ev_annual_demand = {
    "R3_data": 4.4753,
    "R1_data": 5.5913,
    "R2_data": 0.8451,
    "R4_data": 2.4732,
}

vehicle_number = {
    "R3_data": 1714.436,
    "R1_data": 2141.963,
    "R2_data": 323.748,
    "R4_data": 947.455,
}
# %% [markdown]
# Then we deifne the availability of V2G, i.e. the share of electric cars equipped with uni or bi-directional charging
# [%/100] and calculate the respective number of vehicles per model node, which will then be used to set the
# the converter upper limit for the profiles
# Note that all cars without uni- or bi-directional charging control are charged according to the uncontrolled profile.
# %%
smartcharging_share = {"R3_data": 0.8, "R1_data": 0.8, "R2_data": 0.8, "R4_data": 0.8}

ev_uncontrolled_vehicle_number = {
    key: round(value * (1 - smartcharging_share[key]), 4)
    for key, value in vehicle_number.items()
}

ev_controlled_vehicle_number = {
    key: round(value * (smartcharging_share[key]), 4)
    for key, value in vehicle_number.items()
}

ev_parameters = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [m.set.nodesdata, m.set.yearssel, ["EVs_UC", "EVs_CC"]]
    )
)

# number of uncontrolled charged vehicles
ev_parameters.loc[
    idx["R3_data", :, "EVs_UC"], "EVs_Number"
] = ev_uncontrolled_vehicle_number["R3_data"]
ev_parameters.loc[
    idx["R1_data", :, "EVs_UC"], "EVs_Number"
] = ev_uncontrolled_vehicle_number["R1_data"]
ev_parameters.loc[
    idx["R2_data", :, "EVs_UC"], "EVs_Number"
] = ev_uncontrolled_vehicle_number["R2_data"]
ev_parameters.loc[
    idx["R4_data", :, "EVs_UC"], "EVs_Number"
] = ev_uncontrolled_vehicle_number["R4_data"]

# number of controlled charged vehicles
ev_parameters.loc[
    idx["R3_data", :, "EVs_CC"], "EVs_Number"
] = ev_controlled_vehicle_number["R3_data"]
ev_parameters.loc[
    idx["R1_data", :, "EVs_CC"], "EVs_Number"
] = ev_controlled_vehicle_number["R1_data"]
ev_parameters.loc[
    idx["R2_data", :, "EVs_CC"], "EVs_Number"
] = ev_controlled_vehicle_number["R2_data"]
ev_parameters.loc[
    idx["R4_data", :, "EVs_CC"], "EVs_Number"
] = ev_controlled_vehicle_number["R4_data"]

ev_parameters.index.names = ["region", "year", "tech"]
ev_parameters
# %% [markdown]
# EV Storage Converter (Charging Infrastucture - Energy System) for controlled charging (CC)
# First we define all converters needed for CC (controlled charging)

# %%
# technology parameters
converter_techParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([["EVs_CC"], m.set.yearssel])
)
converter_techParam.loc[idx["EVs_CC", :], "lifeTime"] = 15
converter_techParam.loc[idx["EVs_CC", :], "activityUpperLimit"] = 1

m.parameter.add(converter_techParam, "converter_techparam")
converter_techParam
# %%
# converter coefficients
converter_coefficient = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["EVs_CC"], m.set.yearssel, ["Charge", "Discharge"], ["Elec", "EV_Stored"]]
    )
)
converter_coefficient.loc[idx["EVs_CC", :, "Charge", "Elec"], "coefficient"] = -1
converter_coefficient.loc[idx["EVs_CC", :, "Charge", "EV_Stored"], "coefficient"] = 0.99
converter_coefficient.loc[idx["EVs_CC", :, "Discharge", "Elec"], "coefficient"] = 1
converter_coefficient.loc[
    idx["EVs_CC", :, "Discharge", "EV_Stored"], "coefficient"
] = -0.99

converter_coefficient.dropna(how="all", inplace=True)
m.parameter.add(converter_coefficient, "converter_coefficient")
converter_coefficient

# %%
# accounting converter units
accounting_converterUnits = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["Invest", "OMFix"], ["global"], ["EVs_CC"], m.set.yearssel]
    )
).sort_index()
accounting_converterUnits.loc[
    idx["Invest", :, "EVs_CC", :], "perUnitBuild"
] = 0.0  # Mio EUR / unit
accounting_converterUnits.loc[idx["Invest", :, "EVs_CC", :], "useAnnuity"] = 1
accounting_converterUnits.loc[idx["Invest", :, "EVs_CC", :], "amorTime"] = 30
accounting_converterUnits.loc[idx["Invest", :, "EVs_CC", :], "interest"] = 0.06
accounting_converterUnits.loc[idx["OMFix", :, "EVs_CC", :], "perUnitTotal"] = 16.0
accounting_converterUnits.fillna(0, inplace=True)

m.parameter.add(accounting_converterUnits, "accounting_converterunits")
accounting_converterUnits
# %%
# Adding operational costs for the activities of the converter
# (charging infrastructure activities towards energy system)
accounting_converterActivity = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["OMVar"], ["global"], ["EVs_CC"], m.set.yearssel, ["Charge", "Discharge"]]
    )
)
accounting_converterActivity["perActivity"] = 0.0
# Add values (KEur/MWh)
accounting_converterActivity.loc[
    idx["OMVar", :, "EVs_CC", :, "Charge"], "perActivity"
] = 0.001
accounting_converterActivity.loc[
    idx["OMVar", :, "EVs_CC", :, "Discharge"], "perActivity"
] = 0.001

m.parameter.add(accounting_converterActivity, "accounting_converteractivity")
accounting_converterActivity
# %%
# read-in profiles for converters (region, tech, year, type, profile)
# .set_index("cet_cest_timestamp")
charging_avail = pd.read_csv("../_input/vencopy_profiles/vencopy_ChargeAvail_2050.csv")

charging_avail.rename(
    columns={"Unnamed: 0": "tech", "Unnamed: 1": "timestamp"}, inplace=True
)
charging_avail["type"] = "upper"
charging_avail["year"] = "2030"
charging_avail["tech"] = charging_avail["tech"].str.replace("BEV", "EVs_CC")
charging_avail = charging_avail[
    [
        "DE_BadenWue",
        "DE_Bayern",
        "DE_Thueringen",
        "DE_Hessen",
        "type",
        "year",
        "tech",
        "timestamp",
    ]
]
charging_avail.rename(
    columns={
        "DE_BadenWue": "R3_data",
        "DE_Bayern": "R1_data",
        "DE_Thueringen": "R2_data",
        "DE_Hessen": "R4_data",
    },
    inplace=True,
)
charging_avail = charging_avail[charging_avail["tech"] == "EVs_CC"]
charging_avail = charging_avail.melt(id_vars=["tech", "timestamp", "year", "type"])
charging_avail.rename(
    columns={"value": "chargingAvailabililty", "variable": "region"}, inplace=True
)
charging_avail = charging_avail.set_index(
    ["region", "tech", "year", "type", "timestamp"]
).unstack(4)
charging_avail.columns = charging_avail.columns.droplevel(0)

charging_avail.iloc[:, 0:8]
# %%
# scaling the charging availability based on the number of vehicles
charging_numpy = charging_avail.to_numpy()
ev_parameters_numpy = ev_parameters.loc[:, :, "EVs_CC"].to_numpy()

charging_avail_numpy = charging_numpy * ev_parameters_numpy
charging_avail_scaled = pd.DataFrame(
    charging_avail_numpy, columns=charging_avail.columns, index=charging_avail.index
)

charging_avail_scaled_stack = pd.DataFrame(charging_avail_scaled.stack())
max_capacity = charging_avail_scaled_stack.groupby(["region", "year", "tech"]).agg("max")
max_capacity = max_capacity.rename(columns={0: "unitsUpperLimit"})
max_capacity = max_capacity / 1000
max_capacity

# set maximum charging capacity (unitsUpperLimit) based on vehicle number
converter_capacityParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([m.set.nodesdata, m.set.yearssel, ["EVs_CC"]])
)
converter_capacityParam.index.names = ["region", "year", "tech"]
converter_capacityParam["unitsUpperLimit"] = max_capacity["unitsUpperLimit"]
converter_capacityParam["unitsLowerLimit"] = max_capacity["unitsUpperLimit"]

m.parameter.add(converter_capacityParam, "converter_capacityparam")
converter_capacityParam

# %%
# setting the converter activity profile
charging_avail_scaled_numpy = charging_avail_scaled.to_numpy()
max_capacity_numpy = max_capacity.to_numpy()

charging_avail_scaled_numpy = charging_avail_scaled_numpy / max_capacity_numpy
converter_activityProfile = pd.DataFrame(
    charging_avail_scaled_numpy,
    columns=charging_avail_scaled.columns,
    index=charging_avail_scaled.index,
)
converter_activityProfile = converter_activityProfile / 1000
converter_activityProfile = converter_activityProfile.swaplevel(i=1, j=2, axis=0)

m.profile.add(converter_activityProfile, "converter_activityprofile")
converter_activityProfile.iloc[:, 0:8]
# %% [markdown]
# ### Adding a storage reservoir for the battery fleet
#
# Now that we have a converter representing the charging station we need to add a storage reservoir that represent the
# fleet battery
#
# %%
# EV Storage Reservoir
# technology parameters
storage_techParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([["EV_Storage"], m.set.yearssel])
)
storage_techParam.loc[idx["EV_Storage", :], "lifeTime"] = 15
storage_techParam.loc[
    idx["EV_Storage", :], "levelUpperLimit"
] = 1  # does not come into play
storage_techParam.loc[
    idx["EV_Storage", :], "levelLowerLimit"
] = 0  # does not come into play

m.parameter.add(storage_techParam, "storage_techparam")
storage_techParam
# %%
# Annual capacity of EV fleet per node
# R3 4.473 TWh annual = 4473 GWh annual
# R1 5.59 TWh annual = 5590 GWh annual
# R2 2.47 TWh annual = 2470 GWh annual
# R4 0.8451 TWh annual = 845.1 GWh annual

# storage reservoir parameters
storage_reservoirParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([m.set.nodesdata, m.set.yearssel, ["EV_Storage"]])
)
# GWh_el maximum fleet capacity
storage_reservoirParam.loc[idx["R3_data", :, "EV_Storage"], "unitsUpperLimit"] = (
    4473 / 1000
)
storage_reservoirParam.loc[idx["R1_data", :, "EV_Storage"], "unitsUpperLimit"] = (
    5590 / 1000
)
# GWh_el maximum fleet capacity
storage_reservoirParam.loc[idx["R2_data", :, "EV_Storage"], "unitsUpperLimit"] = (
    2470 / 1000
)
storage_reservoirParam.loc[idx["R4_data", :, "EV_Storage"], "unitsUpperLimit"] = (
    845.1 / 1000
)

storage_reservoirParam.loc[idx[:, :, "EV_Storage"], "unitsLowerLimit"] = 0

m.parameter.add(storage_reservoirParam, "storage_reservoirparam")
storage_reservoirParam
# %%
# storage size parameter
storage_sizeParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([["EV_Storage"], m.set.yearssel, ["EV_Stored"]])
)
# GWh_ch (maximum capacity model can use if 1 * unitsUpperLimit see above)
storage_sizeParam.loc[idx["EV_Storage", :, "EV_Stored"], "size"] = 1
storage_sizeParam = storage_sizeParam.dropna()

m.parameter.add(storage_sizeParam, "storage_sizeparam")
storage_sizeParam
# %%
# accounting for the storage unit
accounting_storageUnits = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["Invest", "OMFix"], ["global"], ["EV_Storage"], m.set.yearssel]
    )
)
# Since our storage unit can store 1 GWh we need to scale the Mio  EUR / GWh value
accounting_storageUnits.loc[idx["Invest", :, :, :], "perUnitBuild"] = 0 * 1
accounting_storageUnits.loc[idx["Invest", :, :, :], "useAnnuity"] = 1
accounting_storageUnits.loc[idx["Invest", :, :, :], "amorTime"] = 15
accounting_storageUnits.loc[idx["Invest", :, :, :], "interest"] = 0.06
accounting_storageUnits.loc[idx["OMFix", :, :, :], "perUnitTotal"] = 0 * 1 * 0.015
accounting_storageUnits = accounting_storageUnits.dropna(how="all").fillna(0)

m.parameter.add(accounting_storageUnits, "accounting_storageunits")
accounting_storageUnits
# %%
# ### Adding EV Storage (SoC) profiles
#
# EV Flexible Load
# The flexible load profile describes the fraction of the capacity that can be shifted at the given time,
# it is associated with the converter component of EV storage.

# read-in profiles for storage upper and lower bound (SoC min and max)
SoCMin = pd.read_csv("../_input/vencopy_profiles/vencopy_BatMin_2050.csv")

SoCMin.rename(columns={"Unnamed: 0": "tech", "Unnamed: 1": "timestamp"}, inplace=True)
SoCMin["type"] = "lower"
SoCMin["years"] = "2030"
SoCMin["tech"] = SoCMin["tech"].str.replace("BEV", "EV_Storage")
SoCMin = SoCMin[
    [
        "DE_BadenWue",
        "DE_Bayern",
        "DE_Hessen",
        "DE_Thueringen",
        "type",
        "years",
        "tech",
        "timestamp",
    ]
]
SoCMin.rename(
    columns={
        "DE_BadenWue": "R3_data",
        "DE_Bayern": "R1_data",
        "DE_Thueringen": "R2_data",
        "DE_Hessen": "R4_data",
    },
    inplace=True,
)
SoCMin = SoCMin[SoCMin["tech"] == "EV_Storage"]
SoCMin = SoCMin.melt(id_vars=["tech", "timestamp", "years", "type"])
SoCMin.rename(columns={"value": "SoCMin", "variable": "region"}, inplace=True)
SoCMin = SoCMin.set_index(["region", "years", "tech", "type", "timestamp"]).unstack(4)
SoCMin.columns = SoCMin.columns.droplevel(0)

SoCMax = pd.read_csv("../_input/vencopy_profiles/vencopy_BatMax_2050.csv")
SoCMax.rename(columns={"Unnamed: 0": "tech", "Unnamed: 1": "timestamp"}, inplace=True)
SoCMax["type"] = "upper"
SoCMax["years"] = "2030"
SoCMax["tech"] = SoCMax["tech"].str.replace("BEV", "EV_Storage")
SoCMax = SoCMax[
    [
        "DE_BadenWue",
        "DE_Bayern",
        "DE_Hessen",
        "DE_Thueringen",
        "type",
        "years",
        "tech",
        "timestamp",
    ]
]
SoCMax.rename(
    columns={
        "DE_BadenWue": "R3_data",
        "DE_Bayern": "R1_data",
        "DE_Thueringen": "R2_data",
        "DE_Hessen": "R4_data",
    },
    inplace=True,
)
SoCMax = SoCMax[SoCMax["tech"] == "EV_Storage"]
SoCMax = SoCMax.melt(id_vars=["tech", "timestamp", "years", "type"])
SoCMax.rename(columns={"value": "SoCMax", "variable": "region"}, inplace=True)
SoCMax = SoCMax.set_index(["region", "years", "tech", "type", "timestamp"]).unstack(4)
SoCMax.columns = SoCMax.columns.droplevel(0)
# %%
SoCMin_numpy = SoCMin.to_numpy()
vehicle_number_numpy = np.array(list(vehicle_number.values())).reshape(4, 1)

SoCMin_scaled_numpy = SoCMin_numpy * vehicle_number_numpy
SoCMin_scaled = pd.DataFrame(
    SoCMin_scaled_numpy, columns=SoCMin.columns, index=SoCMin.index
)
SoCMin_scaled_stack = pd.DataFrame(SoCMin_scaled.stack())

SoCMax_numpy = SoCMax.to_numpy()
SoCMax_scaled_numpy = SoCMax_numpy * vehicle_number_numpy
SoCMax_scaled = pd.DataFrame(
    SoCMax_scaled_numpy, columns=SoCMax.columns, index=SoCMax.index
)
SoCMax_scaled_stack = pd.DataFrame(SoCMax_scaled.stack())
max_SoCMax_scaled = SoCMax_scaled_stack.groupby(["region", "years", "tech"]).agg("max")
max_SoCMax_scaled = max_SoCMax_scaled.rename(columns={0: "unitsUpperLimit"})
max_SoCMax_scaled = max_SoCMax_scaled / 1000

SoCMin_scaled

# Rescaling of both min and max profiles
SoCMin_scaled_numpy = SoCMin_scaled.to_numpy()
SoCMax_scaled_numpy = SoCMax_scaled.to_numpy()

max_SoCMax_scaled_numpy = max_SoCMax_scaled.to_numpy()

SoCMax_scaled_numpy = SoCMax_scaled_numpy / max_SoCMax_scaled_numpy
SoCMin_scaled_numpy = SoCMin_scaled_numpy / max_SoCMax_scaled_numpy

SoCMin_minmax = pd.DataFrame(
    SoCMin_scaled_numpy, columns=SoCMin.columns, index=SoCMin.index
)
SoCMax_minmax = pd.DataFrame(
    SoCMax_scaled_numpy, columns=SoCMax.columns, index=SoCMax.index
)

SoCMin_minmax = SoCMin_minmax / 1000
SoCMax_minmax = SoCMax_minmax / 1000
SoCMin_minmax.index.names = [None for i in range(len(SoCMin_minmax.index.names))]
m.profile.add(SoCMin_minmax, "storage_levelprofile")
m.profile.add(SoCMax_minmax, "storage_levelprofile")
SoCMax_minmax.iloc[:, 0:8]
# %% [markdown]
# ### Adding EV Demand for Driving (Sink)

# We now need to add a demand for driving, which represent the power needed for propulsion. It is modelled as a sink
# which is supplied by a new commodity called 'electricity for driving', which represents the outflow from the vehicle
# battery to the electric motor.

# %%
# Converter for driving (converter within the car)
# Define converter, define activity for converter, define indicator for converter and for activity
# Defining a converter, which converts EV stored electricity into electricity for driving (actual demand - fixed)

converter_techParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([["EVs"], m.set.yearssel])
)
converter_techParam.loc[idx[["EVs"], :], "lifeTime"] = 15
converter_techParam.loc[idx[["EVs"], :], "activityUpperLimit"] = 1

m.parameter.add(converter_techParam, "converter_techparam")
converter_techParam
# %%
# converter coefficients
converter_coefficient = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [
            ["EVs"],
            m.set.yearssel,
            ["Discharge_Driving"],
            ["EV_Stored", "driving_energy"],
        ]
    )
)
converter_coefficient.loc[
    idx["EVs", :, "Discharge_Driving", "EV_Stored"], "coefficient"
] = -1.05  # EV_storage to EVs
converter_coefficient.loc[
    idx["EVs", :, "Discharge_Driving", "driving_energy"], "coefficient"
] = 1  # Actual Demand for driving
converter_coefficient.dropna(how="all", inplace=True)

m.parameter.add(converter_coefficient, "converter_coefficient")
converter_coefficient
# %%
# read in profiles for driving demand  (region, tech, year, type, profile)
# .set_index("cet_cest_timestamp")
demand_driving = pd.read_csv("../_input/vencopy_profiles/vencopy_DrivPower_2050.csv")

demand_driving.rename(
    columns={"Unnamed: 0": "tech", "Unnamed: 1": "timestamp"}, inplace=True
)
demand_driving["type"] = "fixed"
demand_driving["years"] = "2030"
demand_driving["commodity"] = "driving_energy"
demand_driving["tech"] = demand_driving["tech"].str.replace("BEV", "DemandEV")
demand_driving = demand_driving[
    [
        "DE_BadenWue",
        "DE_Bayern",
        "DE_Hessen",
        "DE_Thueringen",
        "type",
        "years",
        "tech",
        "timestamp",
        "commodity",
    ]
]
demand_driving = demand_driving[demand_driving["tech"] == "DemandEV"]
demand_driving.rename(
    columns={
        "DE_BadenWue": "R3_data",
        "DE_Bayern": "R1_data",
        "DE_Thueringen": "R2_data",
        "DE_Hessen": "R4_data",
    },
    inplace=True,
)
demand_driving = demand_driving.melt(
    id_vars=["tech", "timestamp", "years", "type", "commodity"]
)
demand_driving.rename(
    columns={"value": "driving_energy", "variable": "region"}, inplace=True
)
demand_driving = demand_driving.set_index(
    ["region", "years", "tech", "commodity", "type", "timestamp"]
).unstack(5)
demand_driving.columns = demand_driving.columns.droplevel(0)

demand_driving

demand_driving_numpy = demand_driving.to_numpy()
vehicle_number_numpy = np.array(list(vehicle_number.values())).reshape(4, 1)

demand_driving_scaled_numpy = -(demand_driving_numpy * vehicle_number_numpy)
demand_driving_scaled = pd.DataFrame(
    demand_driving_scaled_numpy,
    columns=demand_driving.columns,
    index=demand_driving.index,
)
sourcesink_profile = demand_driving_scaled / 1000

m.profile.add(sourcesink_profile, "sourcesink_profile")
sourcesink_profile.iloc[:, 0:5]
# %%
# demand configuration (sets fixed profile for demand for propulsion)
idx_array = [
    list(set(demand_driving.index.get_level_values(i)))
    for i in range(len(demand_driving.index.levels) - 1)
]
df_driving_demand_cfg = pd.DataFrame(index=pd.MultiIndex.from_product(idx_array))
df_driving_demand_cfg.loc[idx[:, :, :, :], "usesFixedProfile"] = 1

m.parameter.add(df_driving_demand_cfg, "sourcesink_config")
df_driving_demand_cfg

# %%
max_capacity = pd.DataFrame({"unitsUpperLimit": -demand_driving_scaled.T.min()})
max_capacity = max_capacity.droplevel(-1).droplevel(-1)
max_capacity.index.names = ["region", "year", "tech"]

max_capacity
# %%
# "converter_capacityParam"
# set max charging capacity based on vehicle number
converter_capacityParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([m.set.nodesdata, m.set.yearssel, ["EVs"]])
)
converter_capacityParam.index.names = ["region", "year", "tech"]
converter_capacityParam["unitsUpperLimit"] = max_capacity["unitsUpperLimit"].values
converter_capacityParam["unitsLowerLimit"] = max_capacity["unitsUpperLimit"].values

m.parameter.add(converter_capacityParam, "converter_capacityparam")
converter_capacityParam
# %% [markdown]
# ### EV Storage Converter (Charging Infrastucture - Energy System) for uncontrolled charging (UC)
#
# Uncontrolled charging is simply modelled as an additional inflexible demand. To this scope, we define a new converter
# (the charging station) which converts the electricity coming from the electrical grid into Li-Ion stored in the
# EV fleet battery.

# %%
# "converter_techParam"
# define UC (uncontrolled charging) converter
# technology parameters
converter_techParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([["EVs_UC"], m.set.yearssel])
)
converter_techParam.loc[idx["EVs_UC", :], "lifeTime"] = 15
converter_techParam.loc[idx["EVs_UC", :], "activityUpperLimit"] = 1

m.parameter.add(converter_techParam, "converter_techparam")
converter_techParam
# %%
# "converter_coefficient"
# converter from energy system to the storage
converter_coefficient = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["EVs_UC"], m.set.yearssel, ["Charge"], ["Elec", "EV_Stored"]]
    )
)
converter_coefficient.loc[idx["EVs_UC", :, "Charge", "Elec"], "coefficient"] = -1
converter_coefficient.loc[
    idx["EVs_UC", :, "Charge", "EV_Stored"], "coefficient"
] = 0.95  # efficiency from grid to battery (eff charging station)
converter_coefficient.dropna(how="all", inplace=True)

m.parameter.add(converter_coefficient, "converter_coefficient")
converter_coefficient
# %%
accounting_converterUnits = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["Invest", "OMFix"], ["global"], ["EVs_UC"], m.set.yearssel]
    )
).sort_index()
accounting_converterUnits.loc[
    idx["Invest", :, "EVs_UC", :], "perUnitBuild"
] = 0.0  # Mio EUR / unit
accounting_converterUnits.loc[idx["Invest", :, "EVs_UC", :], "useAnnuity"] = 1
accounting_converterUnits.loc[idx["Invest", :, "EVs_UC", :], "amorTime"] = 30
accounting_converterUnits.loc[idx["Invest", :, "EVs_UC", :], "interest"] = 0.06
accounting_converterUnits.loc[idx["OMFix", :, "EVs_UC", :], "perUnitTotal"] = 16.0
accounting_converterUnits.fillna(0, inplace=True)

m.parameter.add(accounting_converterUnits, "accounting_converterunits")
accounting_converterUnits
# %%
# "accounting_converterActivity"
# Adding operational costs for the activities of the converter (charging infrastructure activities towards energy
# system)
accounting_converterActivity = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["OMVar"], ["global"], ["EVs_UC"], m.set.yearssel, ["Charge", "Discharge"]]
    )
)
accounting_converterActivity["perActivity"] = 0.0

# Add values (KEur/MWh)
accounting_converterActivity.loc[
    idx["OMVar", :, "EVs_UC", :, "Charge"], "perActivity"
] = 0.001  # kEur/MWh
accounting_converterActivity.loc[
    idx["OMVar", :, "EVs_UC", :, "Discharge"], "perActivity"
] = 0.001  # kEur/MWh

m.parameter.add(accounting_converterActivity, "accounting_converteractivity")
accounting_converterActivity
# %%
# read-in profiles for converters (region, tech, year, type, profile)
charging_avail_uc = pd.read_csv(
    "../_input/vencopy_profiles/vencopy_UncontrCharge_2050.csv"
)

charging_avail_uc.rename(
    columns={"Unnamed: 0": "tech", "Unnamed: 1": "timestamp"}, inplace=True
)
charging_avail_uc["type"] = "fixed"
charging_avail_uc["years"] = "2030"
charging_avail_uc["tech"] = charging_avail_uc["tech"].str.replace("BEV", "EVs_UC")
charging_avail_uc = charging_avail_uc[charging_avail_uc["tech"] == "EVs_UC"]
charging_avail_uc = charging_avail_uc[
    [
        "DE_BadenWue",
        "DE_Bayern",
        "DE_Hessen",
        "DE_Thueringen",
        "type",
        "years",
        "tech",
        "timestamp",
    ]
]
charging_avail_uc.rename(
    columns={
        "DE_BadenWue": "R3_data",
        "DE_Bayern": "R1_data",
        "DE_Thueringen": "R2_data",
        "DE_Hessen": "R4_data",
    },
    inplace=True,
)
charging_avail_uc = charging_avail_uc.melt(
    id_vars=["tech", "timestamp", "years", "type"]
)
charging_avail_uc.rename(
    columns={"value": "chargingAvailabililty", "variable": "region"}, inplace=True
)
charging_avail_uc = charging_avail_uc.set_index(
    ["region", "years", "tech", "type", "timestamp"]
).unstack(4)
charging_avail_uc.columns = charging_avail_uc.columns.droplevel(0)

m.profile.add(charging_avail_uc, "converter_activityprofile")
charging_avail_uc.iloc[:, 0:8]
# %%
# rescale normalized profiles by total number of vehicles
charging_uc_numpy = charging_avail_uc.to_numpy()
ev_parameters_numpy = ev_parameters.loc[:, :, "EVs_UC"].to_numpy()

charging_avail_uc_numpy = charging_uc_numpy * ev_parameters_numpy
charging_avail_uc_scaled = pd.DataFrame(
    charging_avail_uc_numpy,
    columns=charging_avail_uc.columns,
    index=charging_avail_uc.index,
)
charging_avail_uc_scaled_stack = pd.DataFrame(charging_avail_uc_scaled.stack())
max_capacity_uc = charging_avail_uc_scaled_stack.groupby(
    ["region", "years", "tech"]
).agg("max")
max_capacity_uc = max_capacity_uc.rename(columns={0: "unitsUpperLimit"})

max_capacity_uc
# %%
# "converter_capacityParam"
# set max charging capacity based on vehicle number
converter_capacityParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([m.set.nodesdata, m.set.yearssel, ["EVs_UC"]])
)
converter_capacityParam.index.names = ["region", "years", "tech"]
converter_capacityParam["unitsUpperLimit"] = max_capacity_uc["unitsUpperLimit"]
converter_capacityParam["unitsLowerLimit"] = max_capacity_uc["unitsUpperLimit"]
converter_capacityParam = converter_capacityParam / 1000

m.parameter.add(converter_capacityParam, "converter_capacityparam")
converter_capacityParam
# %%
# scale profiles
charging_avail_uc_scaled_numpy = charging_avail_uc_scaled.to_numpy()
max_capacity_uc_numpy = max_capacity_uc.to_numpy()

charging_avail_uc_scaled_numpy = charging_avail_uc_scaled_numpy / max_capacity_uc_numpy
charging_avail_uc_minmax = pd.DataFrame(
    charging_avail_uc_scaled_numpy,
    columns=charging_avail_uc_scaled.columns,
    index=charging_avail_uc_scaled.index,
)

m.profile.add(charging_avail_uc_minmax, "converter_activityprofile")
charging_avail_uc_minmax.iloc[:, 0:9]
# %%
# write all files to `data/` directory
m.write(fileformat="dat")
# %% [markdown]
#
# That's it. We have successfully added electric vehicles to our model,
# specifying a share of the fleet whose charging can be controlled (in contrast
# to uncontrolled charging).
# We can now move to part b of the tutorial and run the optimization.
