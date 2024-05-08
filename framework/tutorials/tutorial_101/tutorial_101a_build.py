# %% [markdown]
# (tutorial_101_label)=
#
# # Tutorial 101 - Converters, sources and sinks
#
# <div style="text-align: center;">
#
# ![Model overview for tutorial 101](../../img/REMix_tutorial101.svg "Model overview for tutorial 101")
#
# Model overview of tutorial 101
#
# </div>
#
# ## Part a: setting up the model
#
# This is the first tutorial to introduce a way to set up a model in REMix. It presents a basic model with four regions
# including renewable energy sources, conventional power plant technologies, an electrical demand and accounting for
# carbon emissions.
#
# For the general structure of REMix tutorials have a look at the README.
#
# We build a first base model to be used in later tutorials to build up on and include other energy system
# components (like energy storage and transfer) as well as technologies (e.g. electric vehicles) and concepts
# (e.g. demand response).

# %% [markdown]
# ### Setting up Python
#
# In this first section, we are importing the Python packages needed to run the model and later exemplary evaluation.
# There are also directories defined where the model data and optimization results will be stored.

# %%
# importing dependencies
import numpy as np
import pandas as pd

from remix.framework.api.instance import Instance

# define often-used shortcut
idx = pd.IndexSlice
# %% [markdown]
# ### General introduction to building models in REMix
#
# For the setup of a model in REMix, preprocessing of data is necessary.
# To do that, the tutorials make use of Pandas DataFrames.
# These are separately set up and collected in lists, before these are being
# written to files that are used as input to the solver.
#
# For the creation of Pandas (pd) DataFrames, we will typically use the
# pd.DataFrame class.
# In addition, we use the pd.MultiIndex.from_product() method to generate a
# multi-index (e.g. three index layers with the first describing the indicator,
# the second describing the indicator used to derive the first indicator and
# third the years).
#
# In the following section, the lists to collect the Pandas DataFrames in are
# initialized in the Instance `m` (as in "model").
# This object is a container in which we will collect all necessary model data.
#
# Not all of the lists initialized with the Instance `m` will be filled in this
# first tutorial.
# This is especially true for storage technologies and energy transfer.
# These two concepts (and more) will be introduced in later tutorials.
#
# One more note: if you do not provide a feature (i.e. fill an empty list),
# REMix will run anyway without that feature but with the other available
# files/features, unless that feature is strictly necessary, like a regional
# mapping.
#
# If you are not yet familiar with the basic functions of Pandas, you can check
# out the 10-minute tutorial in the Pandas documentation:
# https://pandas.pydata.org/pandas-docs/stable/user_guide/10min.html

# %%
# initialize model structure of REMix
m = Instance()

# setting the directory the model data should be written to
# a folder "./data" in the project directory is the default in REMix
m.datadir = "./data"
# %% [markdown]
#
# When printing `m`, you will see all REMix features it includes.
#
# For the purpose of the REMix tutorials, we have prepared some dummy data with
# time profiles that are loaded here.
# %%
# load input data
profiles = pd.read_csv("../_input/profiles.csv", index_col=0)
# %% [markdown]
# ### Defining the model scope
#
# Here is where the model building starts. First of all, we define the model scope.
#
# The model scope describes the fundamental dimensions of the model, e.g. which
# distinct regions and years are modeled.
#
# #### Spatial scope
#
# - `set.nodesdata` : describes the regions for which input data is provided
# such as profiles and capacities for
# power plants
# - `set.nodesmodel` : describes the model regions which can be the same as the
# data regions if the optimization should be done in full resolution.
# - `map.aggregatenodesmodel` : describes the aggregation mapping for data to
# model regions. This can be a 1:1 mapping (like `R3_data` to `R3_model`) or a
# n:1 mapping (like e.g. "R1_North_data" and "R1_South_data" to `R1_model`) if
# multiple data regions should be summed up to a model region.
#
# #### Temporal scope
#
# - `set.years` : the individual years which can be modeled for historical and
# new power plants
# - `set.yearssel` : the years which should be optimized during the run. For
# now, we only use a single year to be optimized.
#
# Our model will comprise four regions, also referred to as "nodes", whose names
# can be arbitrarily chosen. Here, they are called `R3_model`, `R1_model`,
# `R2_model` and `R4_model` (although having nothing to do with the actual
# energy systems of the countries these abbreviations hint at).
# In the first two tutorials, we will only use one node, which is
# `R1_data`/`R1_model`, so the other nodes are not needed until tutorial 103.

# %%
# "map_aggregateNodesModel"
# DataFrame for aggregation from data to model regions
df = pd.DataFrame(
    [
        ["R1_data", "R1_model", 1],
        ["R2_data", "R2_model", 1],  # not strictly necessary for tutorial 1 and 2
        ["R3_data", "R3_model", 1],  # not strictly necessary for tutorial 1 and 2
        ["R4_data", "R4_model", 1],  # not strictly necessary for tutorial 1 and 2
    ]
)
df.columns = ["nodesData", "nodesModel", "aggregate"]
df.set_index(["nodesData", "nodesModel"], inplace=True)
df["aggregate"] = ""
df.columns = [""]

m.map.add(df, "aggregatenodesmodel")

# Get the data and model regions based on the mapping
# "set_nodesData"
m.set.add(
    list(sorted(set(m.map.aggregatenodesmodel.index.get_level_values(0)))), "nodesdata"
)
# "set_nodesModel" & "set_nodesModelSel"
m.set.add(
    list(sorted(set(m.map.aggregatenodesmodel.index.get_level_values(1)))), "nodesmodel"
)

# Set the years to be considered in the model and the years to be optimized
# "set_years"
m.set.add(
    ["2030"], "years"
)  # must include all years that data is provided for in the model
# "set_yearsSel"
m.set.add(["2030"], "yearssel")  # years to be optimised
# %% [markdown]
# ### Setting the objective function and indicator bounds
#
# Models in REMix are usually optimized based on a cost-minimization approach.
# The framework theoretically also allows other approaches.
#
# We will use different types of commodities - electricity, methane, carbon
# dioxide - and system costs as indicator.
# We will use the following units for these:
#
# - Elec : electricity in GWh_el
# - CH4 : methane in GWh_ch
# - CO2 : carbon dioxide emissions in tsd. t or kt
# - Cost (Invest, OMVar, OMFix, CarbonCost, FuelCost) : cost values in million EUR or MEUR
#
# In the first DataFrame we define a value for the indicator `SystemCost` and
# column `obj` to -1 to communicate that we want to minimize this indicator.
# Similarly, a value of 1 would indicate a maximization.
# The first field is used for the regional and year dimensions.
# The value `global` uses all the regions in the system (in this example
# R1_model, R2_model, R3_model, R4_model), whereas the value `horizon` takes
# into account all years in the set `set.yearssel` (here only 2030).
#
# We set a social discount rate in the same DataFrame, which will be the default
# value throughout the model, but can be overwritten for certain technologies or
# model regions if wanted.

# %%
# "accounting_indicatorBounds"
# setting the objective function and indicator bounds
accounting_indicatorBounds = pd.DataFrame(
    index=pd.MultiIndex.from_product([["global"], ["horizon"], ["SystemCost"]])
)
accounting_indicatorBounds["obj"] = -1  # minimization of system costs
accounting_indicatorBounds["discount"] = 0.07  # social discount rate for the indicators

m.parameter.add(accounting_indicatorBounds, "accounting_indicatorbounds")
accounting_indicatorBounds
# %% [markdown]
# We are also setting up the indicators we want to account for as `SystemCost`
# in the model.
#
# Indicators are used for general accounting inside the energy system. For this
# purpose we introduce an indicator `SystemCost` to reflect the overall costs of
# the system.
# This indicator is calculated by summing up the following individual cost
# indicators with an equal weighting of 1 in the `accounting_perIndicator`
# DataFrame.
#
# - `Invest` : investment cost for a technology unit (in MEUR/MW)
# - `OMVar` : variable operation and maintenance cost  (in MEUR/MWh) (not set in this tutorial)
# - `OMFix` : fix operation and maintenance costs (in MEUR/MW/year)
# - `FuelCost` : costs for imports of methane into the model regions (in MEUR/MWh)

# %%
# "accounting_perIndicator"
# set up accounting per indicator for all years to calculate
accounting_perIndicator = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [
            ["SystemCost"],
            [
                "Invest",
                "OMFix",
                "FuelCost",
            ],
            ["global"],
            m.set.yearssel,  # accounting for all optimization years
        ]
    )
)
accounting_perIndicator["perIndicator"] = 1

m.parameter.add(accounting_perIndicator, "accounting_perindicator")
accounting_perIndicator
# %% [markdown]
# ### Converter technologies
#
# #### Adding converter technologies
#
# In this section, the basic structure of including different converter
# technologies in REMix is introduced.
#
# In this basic model, we introduce the possibility for the model to build
# methane-fired combined-cycle gas turbines ("CCGT"), solar power plants ("PV")
# and onshore wind turbine ("WindOnshore").
#
# The names chosen for the technologies are completely arbitrary.
# We are trying to use the same ones throughout the tutorials, however.

# %%
# "converter_techParam"
# setting technology parameters
converter_techParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([["CCGT", "PV", "WindOnshore"], m.set.yearssel])
)
converter_techParam.loc[idx["CCGT"], "lifeTime"] = 30  # years
converter_techParam.loc[
    idx["CCGT"], "activityUpperLimit"
] = 1  # availability of technology

converter_techParam.loc[idx["PV"], "lifeTime"] = 20  # years
converter_techParam.loc[
    idx["PV"], "activityUpperLimit"
] = 0  # this value will be replaced later on with the normalized feed-in profile

converter_techParam.loc[idx["WindOnshore"], "lifeTime"] = 25
converter_techParam.loc[
    idx["WindOnshore"], "activityUpperLimit"
] = 0  # this value will be replaced later on with the normalized feed-in profile

m.parameter.add(converter_techParam, "converter_techparam")
converter_techParam
# %%
# "converter_capacityParam"
# defining upper and/or lower limits for converter technologies
converter_capacityParam = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [m.set.nodesdata, m.set.yearssel, ["CCGT", "PV", "WindOnshore"]]
    )
)
converter_capacityParam.loc[idx["R1_data", :, "CCGT"], "unitsUpperLimit"] = 100  # GW_el
converter_capacityParam.loc[idx["R1_data", :, "CCGT"], "unitsLowerLimit"] = 0  # GW_el
converter_capacityParam.loc[idx["R1_data", :, "PV"], "unitsUpperLimit"] = 40  # GW_el
converter_capacityParam.loc[idx["R1_data", :, "PV"], "unitsLowerLimit"] = 0  # GW_el
converter_capacityParam.loc[
    idx["R1_data", :, "WindOnshore"], "unitsUpperLimit"
] = 60  # GW_el
converter_capacityParam.loc[
    idx["R1_data", :, "WindOnshore"], "unitsLowerLimit"
] = 0  # GW_el
converter_capacityParam.dropna(how="all", inplace=True)

m.parameter.add(converter_capacityParam, "converter_capacityparam")
converter_capacityParam
# %% [markdown]
# Activities in REMix are the conversion processes a technology can perform.
# For this example we define an activity "Powergen" (as in power generation).
#
# For the CCGT technology this means burning methane in order to get electricity
# and carbon dioxide as a by-product of the combustion process.
#
# For the renewable energy sources wind and PV we model the activity `Powergen`
# by setting a value of 1, which is arbitrary in this case, however, since the
# actual potential for wind and solar energy is modeled as "activityProfile"
# below, which overwrites this value.

# %%
# "converter_coefficient"
converter_coefficient = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [
            ["CCGT", "PV", "WindOnshore"],
            m.set.yearssel,
            ["Powergen"],
            ["CH4", "Elec", "CO2"],
        ]
    )
)
converter_coefficient.loc[idx["CCGT", :, :, "Elec"], "coefficient"] = 1  # GWh_el
converter_coefficient.loc[idx["CCGT", :, :, "CH4"], "coefficient"] = -1.587  # GWh_ch
converter_coefficient.loc[idx["CCGT", :, :, "CO2"], "coefficient"] = 0.320  # kt CO2

converter_coefficient.loc[idx["PV", :, :, "Elec"], "coefficient"] = 1  # GWh_el

converter_coefficient.loc[idx["WindOnshore", :, :, "Elec"], "coefficient"] = 1  # GWh_el
converter_coefficient.dropna(how="all", inplace=True)

m.parameter.add(converter_coefficient, "converter_coefficient")
converter_coefficient
# %% [markdown]
# Since we now introduced a conversion unit that runs on variable renewable
# energy, we need to limit the profile for the activity on the potential
# feed-in.
# We can do this in a similar way to adding the electrical demand profile.
#
# The values in the `profiles.csv` are given in mega watt (MW) of electrical
# feed-in.
# We need to normalize them to values between 0 and 1.
# This normalized profile describes the maximum activity per unit of power plant.
#
# Example: 10 PV units with 1 GW rated capacity each (as specified by the
# activity parameter) with an activity profile of 0.24 in hour 11 could produce
# up to 10 * 1 GWh/h * 0.24 = 2.4 GWh/h.

# %%
# "converter_activityProfile"
# load the profiles DataFrame, select its PV and WindOnshore columns
converter_activityProfile = profiles[["PV", "WindOnshore"]]

# convert from MW to GW
converter_activityProfile = converter_activityProfile.div(1e3).T

converter_activityProfile = converter_activityProfile.div(
    converter_activityProfile.max(axis=1), axis=0
)
converter_activityProfile.index.names = ["techs"]

# add columns and set them as index
converter_activityProfile["region"] = "R1_data"
converter_activityProfile["years"] = "2030"
converter_activityProfile["type"] = "upper"
converter_activityProfile = converter_activityProfile.reset_index().set_index(
    ["region", "years", "techs", "type"]
)

m.profile.add(converter_activityProfile, "converter_activityprofile")
converter_activityProfile.iloc[:, 0:8]
# %%
# "accounting_converterUnits"
# setting the costs of technologies
accounting_converterUnits = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [
            ["Invest", "OMFix"],
            ["global"],
            ["CCGT", "PV", "WindOnshore"],
            m.set.yearssel,
        ]
    )
).sort_index()

accounting_converterUnits.loc[
    idx["Invest", "global", "CCGT", "2030"], "perUnitBuild"
] = 700.0  # Mio EUR per unit
accounting_converterUnits.loc[
    idx["Invest", "global", "CCGT", "2030"], "useAnnuity"
] = 1  # binary yes/no
accounting_converterUnits.loc[
    idx["Invest", "global", "CCGT", "2030"], "amorTime"
] = 30  # years
accounting_converterUnits.loc[
    idx["Invest", "global", "CCGT", "2030"], "interest"
] = 0.06  # percent/100
accounting_converterUnits.loc[
    idx["OMFix", "global", "CCGT", "2030"], "perUnitTotal"
] = 28.0  # Mio EUR per unit

accounting_converterUnits.loc[
    idx["Invest", "global", "PV", "2030"], "perUnitBuild"
] = 518.0
accounting_converterUnits.loc[idx["Invest", "global", "PV", "2030"], "useAnnuity"] = 1
accounting_converterUnits.loc[idx["Invest", "global", "PV", "2030"], "amorTime"] = 20
accounting_converterUnits.loc[idx["Invest", "global", "PV", "2030"], "interest"] = 0.06
accounting_converterUnits.loc[
    idx["OMFix", "global", "PV", "2030"], "perUnitTotal"
] = 7.7

accounting_converterUnits.loc[
    idx["Invest", "global", "WindOnshore", "2030"], "perUnitBuild"
] = 1368.0
accounting_converterUnits.loc[
    idx["Invest", "global", "WindOnshore", "2030"], "useAnnuity"
] = 1
accounting_converterUnits.loc[
    idx["Invest", "global", "WindOnshore", "2030"], "amorTime"
] = 25
accounting_converterUnits.loc[
    idx["Invest", "global", "WindOnshore", "2030"], "interest"
] = 0.025
accounting_converterUnits.loc[
    idx["OMFix", "global", "WindOnshore", "2030"], "perUnitTotal"
] = 25.8
accounting_converterUnits.fillna(0, inplace=True)

m.parameter.add(accounting_converterUnits, "accounting_converterunits")
accounting_converterUnits
# %% [markdown]
# ### Sources and sinks
#
# #### Adding a demand profile as sink
#
# In this part, we set a demand for the data node `R1_data` (which is aggregated
# to the model node `R1_model`) only.
# The region name and year have to be included in the `map.aggregatenodesmodel`
# and `set.years` defined in the beginning.
# The name for the source-sink technology (here: `Demand`) can be freely chosen.
#
# We need to specify that the demand is applied to the electrical commodity and
# that this profile needs to be matched exactly on an hour-by-hour level.

# %%
# "sourcesink_profile"
# load the profiles DataFrame, select the demand column
sourcesink_profile = profiles[["demand_R1"]]

# divide by 1000 to convert to GW, multiply with -1 because this is the
# REMix convention for accounting for sinks/demand
sourcesink_profile = sourcesink_profile.div(1e3).mul(-1)
# transpose DataFrame for needed format
sourcesink_profile = sourcesink_profile.T

# add columns and set them as index
sourcesink_profile["nodesData"] = "R1_data"
sourcesink_profile["years"] = "2030"
sourcesink_profile["techs"] = "Demand"
sourcesink_profile["commodity"] = "Elec"
sourcesink_profile["type"] = "fixed"
sourcesink_profile = sourcesink_profile.set_index(
    ["nodesData", "years", "techs", "commodity", "type"]
)

m.profile.add(sourcesink_profile, "sourcesink_profile")
sourcesink_profile.iloc[:, 0:8]
# %% [markdown]
# Now that we have created the profile, we need to create a config with the
# information that the created profile is going to be integrated into the model
# as fixed profile.

# %%
# "sourcesink_config" (demand configuration)
sourcesink_config = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [m.set.nodesdata, m.set.yearssel, ["Demand"], ["Elec"]]
    )
)
sourcesink_config.loc[idx["R1_data", :, :, :], "usesFixedProfile"] = 1
sourcesink_config.dropna(inplace=True)

m.parameter.add(sourcesink_config, "sourcesink_config")
sourcesink_config
# %% [markdown]
# #### Add sources for fuels and sinks for carbon emissions
#
# Since CCGT uses CH4 as a fuel, we need to allow import of CH4 for the model
# region `R1_model` (since the technology is only installed there).
# This is very similar to the source-sink technology we used for the electrical
# demand.
# However, in this case we want to be able to import an unlimited amount of fuel
# at a fixed price of 0.0306 million EUR/GWh_ch.
# By adding a lower profile of 0, we ensure the model cannot export fuel to make
# money.

# %%
# "sourcesink_annualSum"
# limiting the annual sum of fuel imports into a model region
sourcesink_annualSum = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [m.set.nodesdata, m.set.yearssel, ["FuelImport"], ["CH4"]]
    )
)
sourcesink_annualSum.loc[idx["R1_data", :, :, :], "upper"] = np.inf
sourcesink_annualSum.loc[idx["R1_data", :, :, :], "lower"] = 0
sourcesink_annualSum.dropna(inplace=True)

m.parameter.add(sourcesink_annualSum, "sourcesink_annualsum")
sourcesink_annualSum
# %%
# "sourcesink_config" (import configuration)
sourcesink_config = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [m.set.nodesdata, m.set.yearssel, ["FuelImport"], ["CH4"]]
    )
)
sourcesink_config.loc[idx["R1_data", :, :, :], "usesUpperSum"] = 1
sourcesink_config.loc[idx["R1_data", :, :, :], "usesLowerProfile"] = 1
sourcesink_config.dropna(inplace=True)

m.parameter.add(sourcesink_config, "sourcesink_config")
sourcesink_config
# %%
# "accounting_sourcesinkFlow"
# setting a cost for methane imports
accounting_sourcesinkFlow = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["FuelCost"], ["global"], m.set.yearssel, ["FuelImport"], ["CH4"]]
    )
)
accounting_sourcesinkFlow["perFlow"] = 0.03060  # Mio EUR per GWh_ch CH4

m.parameter.add(accounting_sourcesinkFlow, "accounting_sourcesinkflow")
accounting_sourcesinkFlow
# %% [markdown]
# Similar to the fuel source we need to specify a sink for our carbon emissions.
# In this case we need to use negative values since the carbon is leaving our
# frame of accounting. So we specify a lower sum of -infinity and an upper
# profile of 0 (meaning we are not allowed to extract carbon out of the
# atmosphere).
# By changing the condition from -infinity to -100, we could also impose a
# carbon limit of 100 kilotonnes of CO2.
# Or we could add a new indicator "CarbonCost" (at the top) which accounts for
# the carbon flow out of the system and imposes an associated cost.

# %%
# "sourcesink_annualSum"
# limiting annual sum of carbon emissions
sourcesink_annualSum = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [m.set.nodesdata, m.set.yearssel, ["Emission"], ["CO2"]]
    )
)
sourcesink_annualSum.loc[idx["R1_data", :, :, :], "lower"] = -np.inf
sourcesink_annualSum.dropna(inplace=True)

m.parameter.add(sourcesink_annualSum, "sourcesink_annualsum")
sourcesink_annualSum
# %%
# "sourcesink_config" (emission configuration)
sourcesink_config = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [m.set.nodesdata, m.set.yearssel, ["Emission"], ["CO2"]]
    )
)
sourcesink_config.loc[idx["R1_data", :, :, :], "usesLowerSum"] = 1
sourcesink_config.loc[idx["R1_data", :, :, :], "usesUpperProfile"] = 1
sourcesink_config.dropna(inplace=True)

m.parameter.add(sourcesink_config, "sourcesink_config")
sourcesink_config
# %% [markdown]
# ### Writing DataFrames to files
#
# In this section we collect the DataFrames from the previous sections and
# convert them to files inside the folder which was specified in the
# beginning or the `data/` directory by default. Writing to `*.csv` files will
# work similarly.

# %%
# write all files to the datadir
m.write(fileformat="dat")
# %% [markdown]
# We have finished building our REMix data model now. In part b of this
# tutorial, we are looking into how we can execute it.
