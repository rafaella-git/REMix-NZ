# %% [markdown]
# (tutorial_202_label)=
#
# # Tutorial 202 - Path optimization with carbon budget
#
# <div style="text-align: center;">
#
# ![Model overview for tutorial 202](../../img/REMix_tutorial202.svg "Model overview for tutorial 202")
#
# Model overview of tutorial 202
#
# </div>
#
# ## Part a: setting up the model
#
# %% [markdown]
# ### Setting up Python

# %%
# importing dependencies
import numpy as np
import pandas as pd
from remix.framework import Instance
from remix.framework import read_dat

# define often-used shortcut
idx = pd.IndexSlice
# %% [markdown]
# ### Initialize the model in REMix
#
# %%
# initialize model structure of REMix
m = Instance()

m.datadir = "./data"
# %% [markdown]
# ### Defining the model scope
#
# Let's start building our model!
# First of all, we define the model scope.
#
# The model scope describes the fundamental dimensions of the model, e.g. which
# distinct regions and years are modeled.
#
# %%
# "map_aggregateNodesModel"
# DataFrame for aggregation from data to model regions
df = pd.DataFrame(
    [
        ["R1_data", "R1_model", 1],
        ["R2_data", "R2_model", 1],
        ["R3_data", "R3_model", 1],
        ["R4_data", "R4_model", 1],
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
m.set.add(
    list(sorted(set(m.map.aggregatenodesmodel.index.get_level_values(1)))),
    "nodesmodelsel",
)

# Set the years to be considered in the model and the years to be optimized
# "set_years"
m.set.add(
    ["2020", "2030", "2040", "2050"], "years"
)  # must include all years that data is provided for in the model
# "set_yearsSel"
m.set.add(["2020", "2030", "2040", "2050"], "yearssel")  # years to be optimized
# %% [markdown]
# ### Setting the objective function and indicator bounds
#
# Models in REMix are usually optimized based on a cost-minimization approach.
# The framework theoretically also allows other approaches.
#
# We will use different types of commodities - electricity, methane, carbon
# dioxide - and system costs as indicator.
# We use the following units for these:
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
# into account all years in the set `set.yearssel`.
#
# The `discount` column in this dataframe represents the `social discount rate`
# which represents a societal valuation of the future value.
#
# It is a fundamental component of path optimization, as it sets the framework
# to make decisions between present and future value. This discount is different
# from the discount of the WACC in the sense that the investment decisions done
# in a certain year are subject to the `social discount rate` after being
# calculated.

# %%
# "accounting_indicatorBounds"
# setting the objective function and indicator bounds
accounting_indicatorBounds_global = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["global"], ["horizon"], ["SystemCost", "Carbon", "Slack_Carbon"]]
    )
)

delta_indexes = ["Delta_PV", "Delta_Wind"]

slack_indexes = ["Slack_Delta_PV", "Slack_Delta_Wind"]

accounting_indicatorBounds_global.loc[
    idx[:, :, "SystemCost"], "obj"
] = -1  # minimization of system costs
accounting_indicatorBounds_global.loc[
    idx[:, :, "SystemCost"], "endyear"
] = 10  # duration of the last period (used for discounting cost)
accounting_indicatorBounds_global.loc[
    idx[:, :, "SystemCost"], "integral"
] = 1  # if enabled, indicator will also be accounted for years between yearssel
accounting_indicatorBounds_global.loc[
    idx[:, :, "SystemCost"], "discount"
] = 0.01  # social discount rate

# Emissions of Germany, Italy, Switzerland and France in 2019:
# 1388270.0 kilotonnes/a
# Calculated from: https://edgar.jrc.ec.europa.eu/report_2020
# From which 80% of these emissions across 2020, 2030, 2040 and 2050:
# 4442464.0 kilotonnes
accounting_indicatorBounds_global.loc[
    idx[:, :, "Carbon"], "upperValue"
] = 8000000.0  # rounded value
accounting_indicatorBounds_global.loc[
    idx[:, :, "Carbon"], "endyear"
] = 10  # duration of the last period (used for discounting emissions)
accounting_indicatorBounds_global.loc[
    idx[:, :, "Carbon"], "integral"
] = 1  # if enabled, indicator will also be accounted for years between yearssel
accounting_indicatorBounds_global.loc[
    idx[:, :, "Carbon"], "useUpper"
] = 1  # parameter to indicate the model to use the `upperValue` field

accounting_indicatorBounds_global.loc[
    idx[:, :, "Slack_Carbon"], "isVariable"
] = 1  # parameter to indicate that this is a free variable
accounting_indicatorBounds_global.loc[
    idx[:, :, "Slack_Carbon"], "useLower"
] = 1  # use the lower level of this variable

accounting_indicatorBounds_deltas = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [m.set.nodesdata, m.set.years, slack_indexes + delta_indexes]
    )
)
# adjust Slack variables
accounting_indicatorBounds_deltas.loc[idx[:, :, slack_indexes], "isVariable"] = 1
accounting_indicatorBounds_deltas.loc[idx[:, :, slack_indexes], "useLower"] = 1
# adjust Delta variables
accounting_indicatorBounds_deltas.loc[idx[:, :, delta_indexes], "useUpper"] = 1
accounting_indicatorBounds_deltas.loc[
    idx[["R1_data"], :, delta_indexes], "upperValue"
] = 45
accounting_indicatorBounds_deltas.loc[
    idx[["R2_data"], :, delta_indexes], "upperValue"
] = 60
accounting_indicatorBounds_deltas.loc[
    idx[["R3_data"], :, delta_indexes], "upperValue"
] = 50
accounting_indicatorBounds_deltas.loc[
    idx[["R4_data"], :, delta_indexes], "upperValue"
] = 50


accounting_indicatorBounds = pd.concat(
    [accounting_indicatorBounds_global, accounting_indicatorBounds_deltas]
).fillna(0)

m.parameter.add(accounting_indicatorBounds, "accounting_indicatorbounds")
accounting_indicatorBounds
# %% [markdown]
# We are also setting up the indicators we want to account for as `SystemCost`
# in the model.
#
# Indicators are used for general accounting inside the energy system.
# For this purpose we introduce an indicator `SystemCost` to reflect the overall
# costs of the system. This indicator is calculated by summing up the following
# individual cost indicators with an equal weighting of 1 in the
# `accounting_perIndicator` DataFrame.
#
# - `Invest` : investment cost for a technology unit (in MEUR/MW)
# - `OMFix` : fix operation and maintenance costs (in MEUR/MW/year)
# - `OMVar` : variable operation and maintenance cost (in MEUR/MWh)
# - `FuelCost` : costs for imports of methane into the model regions (in MEUR/MWh)
#
# Indicators can also be used to set soft constraints in the system like how
# much of a certain technology can be built within a timespan.
# We account for this by using:
#
# - `Delta_PV` : sets an upper limit of the 10-year installation of new PV
# capacity.
# - `Delta_Wind` : sets an upper limit of the 10-year installation of new
# wind-turbine capacity.
#
# Other use for indicators is to track if the constraints built in the system
# are to tight. To do that, we define slack indicators which allow the model to
# be solvable in such cases:
#
# - `Slack_Carbon` : assigns a cost to excess emissions of `Carbon`. This can be
# understood as carbon tax.
# - `Slack_Delta_PV` : assigns an additional cost if more capacity is needed
# than allowed with `Delta_PV` (prevents infeasible model).
# - `Slack_Delta_Wind` : assigns an additional cost if more capacity is needed
# than allowed with `Delta_Wind` (prevents infeasible model).
#
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
                "OMVar",
                "FuelCost",
                "Slack_Carbon",
                "Slack_Delta_PV",
                "Slack_Delta_Wind",
            ],
            ["global"],
            ["horizon"],  # "horizon" := accounting for all optimization years
        ]
    )
)
accounting_perIndicator.loc[
    idx[["SystemCost"], ["Invest", "OMFix", "OMVar", "FuelCost"]], "perIndicator"
] = 1
accounting_perIndicator.loc[idx[["SystemCost"], slack_indexes], "perIndicator"] = 1500
accounting_perIndicator.loc[idx[["SystemCost"], ["Slack_Carbon"]], "perIndicator"] = 0.6

accounting_perIndicator_slack = pd.DataFrame(
    {"perIndicator": [-1, -1, -1]},
    index=pd.MultiIndex.from_tuples(
        [
            ("Carbon", "Slack_Carbon", "global", "horizon"),
            ("Delta_PV", "Slack_Delta_PV", "global", "horizon"),
            ("Delta_Wind", "Slack_Delta_Wind", "global", "horizon"),
        ]
    ),
)
accounting_perIndicator = pd.concat(
    [accounting_perIndicator, accounting_perIndicator_slack]
)

m.parameter.add(accounting_perIndicator, "accounting_perindicator")
accounting_perIndicator
# %%
# "accounting_converterActivity"
# set up accounting per indicator for all years to calculate
accounting_converteractivity = read_dat("input/accounting_converteractivity.dat")
m.parameter.add(accounting_converteractivity, "accounting_converteractivity")
accounting_converteractivity
# %% [markdown]
# ### Converter technologies
#
# #### Adding converter technologies
#
# In this section, the basic structure of including different converter
# technologies in REMix is introduced.
#
# In this basic model, we introduce the possibility for the model to build
# methane-fired combined-cycle gas turbines (CCGT), regular gas turbines
# (OCGT), solar power plants (PV) and onshore wind turbine (WindOnshore).
# The names chosen for the technologies are completely arbitrary.
#
# There is two special converter technologies in this model. `GreenImport`
# is going to handle the gas imports from Power2Gas and `LiIon` is the charger
# of the lithium-ion storage.
# %%
# "converter_techParam"
# setting technology parameters
m.set.add(
    ["CCGT", "OCGT", "LiIon", "PV", "WindOnshore", "GreenImport"],
    "converter_techs",
)
converter_techParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([m.set.converter_techs, ["2020"]])
)
converter_techParam.loc[
    idx[["CCGT", "OCGT", "GreenImport"], :], "lifeTime"
] = 30  # years
converter_techParam.loc[
    idx[["CCGT", "OCGT", "GreenImport", "LiIon"], :], "activityUpperLimit"
] = 1  # availability of technology

converter_techParam.loc[idx["PV"], "lifeTime"] = 20  # years
converter_techParam.loc[
    idx["PV"], "activityUpperLimit"
] = 0  # this value will be replaced later on with the normalized feed-in profile

converter_techParam.loc[idx[["WindOnshore", "LiIon"], :], "lifeTime"] = 25
converter_techParam.loc[
    idx["WindOnshore"], "activityUpperLimit"
] = 0  # this value will be replaced later on with the normalized feed-in profile

m.parameter.add(converter_techParam, "converter_techparam")
converter_techParam
# %%
# "converter_capacityParam"
# defining upper and/or lower limits for converter technologies
converter_capacityParam = read_dat("input/converter_capacityparam.dat")

m.parameter.add(converter_capacityParam, "converter_capacityparam")
converter_capacityParam
# %% [markdown]
# Activities in REMix are the conversion processes a technology can perform.
# For this example we define the activities "Powergen" (as in power generation),
# "Charge", "Discharge" and "Import".
#
# CCGT and OCGT burn methane in order to get electricity and carbon dioxide as a
# by-product of the combustion process for their "Powergen" activity.
#
# We also need to set an activity for renewable energy sources for the activity
# `Powergen`. This value is arbitrary in our case, however, since it is
# overwritten by their actual potentials below (as "activityProfile", see below).
#
# The `Import` activity will handle the gas supply of GreenImport.
# %%
# "converter_coefficient"
converter_coefficient = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [
            m.set.converter_techs,
            ["2020"],
            ["Powergen", "Charge", "Discharge", "Import"],
            ["CH4", "Elec", "CO2", "Elec_LiIon"],
        ]
    )
)
converter_coefficient.loc[
    idx["CCGT", :, "Powergen", "Elec"], "coefficient"
] = 1  # GWh_el
converter_coefficient.loc[
    idx["CCGT", :, "Powergen", "CH4"], "coefficient"
] = -1.59  # GWh_ch

converter_coefficient.loc[
    idx["OCGT", :, "Powergen", "Elec"], "coefficient"
] = 1  # GWh_el
converter_coefficient.loc[
    idx["OCGT", :, "Powergen", "CH4"], "coefficient"
] = -2.13  # GWh_ch

converter_coefficient.loc[idx["PV", :, "Powergen", "Elec"], "coefficient"] = 1  # GWh_el
converter_coefficient.loc[
    idx["WindOnshore", :, "Powergen", "Elec"], "coefficient"
] = 1  # GWh_el
converter_coefficient.loc[
    idx["GreenImport", :, "Import", "CH4"], "coefficient"
] = 1  # GWh_ch

converter_coefficient.loc[
    idx["LiIon", :, "Charge", "Elec"], "coefficient"
] = -1  # GWh_el
converter_coefficient.loc[
    idx["LiIon", :, "Charge", "Elec_LiIon"], "coefficient"
] = 0.95  # GWh_el

converter_coefficient.loc[
    idx["LiIon", :, "Discharge", "Elec"], "coefficient"
] = 0.95  # GWh_el
converter_coefficient.loc[
    idx["LiIon", :, "Discharge", "Elec_LiIon"], "coefficient"
] = -1  # GWh_el

converter_coefficient.dropna(how="all", inplace=True)

m.parameter.add(converter_coefficient, "converter_coefficient")
converter_coefficient
# %% [markdown]
# Since we introduced a conversion unit that runs on variable renewable energy,
# we need to limit the profile for the activity on the potential feed-in.
#
# The values in the `profiles.csv` are given in mega watt (MW) of electrical
# feed-in. We need to normalize them to values between 0 and 1.
# This normalized profile describes the maximum activity per unit of power plant.
#
# Example: 10 PV units with 1 GW rated capacity each (as specified by the
# activity parameter) with an activity profile of 0.24 in hour 11 could produce
# up to 10 * 1 GWh/h * 0.24 = 2.4 GWh/h.

# %%
# "converter_activityProfile"
# load the profiles DataFrame, select its PV and WindOnshore columns
profiles = pd.read_csv("../_input/profiles.csv", index_col=0)

for data_node in m.set.nodesdata:
    for year in m.set.yearssel:
        converter_activityProfile = profiles[["PV", "WindOnshore"]]

        # convert from MW to GW
        converter_activityProfile = converter_activityProfile.div(1e3).T

        converter_activityProfile = converter_activityProfile.div(
            converter_activityProfile.max(axis=1), axis=0
        )
        converter_activityProfile.index.names = ["techs"]

        converter_activityProfile["region"] = data_node
        converter_activityProfile["years"] = year
        converter_activityProfile["type"] = "upper"
        converter_activityProfile = converter_activityProfile.reset_index().set_index(
            ["region", "years", "techs", "type"]
        )

        m.profile.add(converter_activityProfile, "converter_activityprofile")

m.profile.converter_activityprofile
converter_activityProfile.iloc[:, 0:8]

m.profile.converter_activityprofile.iloc[:, 0:8]
# %%
# "accounting_converterUnits"
# setting the costs of technologies
accounting_converterUnits = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [
            ["Invest", "OMFix", "Delta_PV", "Delta_Wind"],
            ["global"],
            m.set.converter_techs,
            m.set.yearssel,
        ]
    )
).sort_index()

accounting_converterUnits.loc[
    idx["Invest", "global", "CCGT", "2020"], "perUnitBuild"
] = 700.0  # Mio EUR per unit
accounting_converterUnits.loc[
    idx["Invest", "global", "CCGT", "2020"], "useAnnuity"
] = 1  # binary yes/no
accounting_converterUnits.loc[
    idx["Invest", "global", "CCGT", "2020"], "amorTime"
] = 30  # years
accounting_converterUnits.loc[
    idx["Invest", "global", "CCGT", "2020"], "interest"
] = 0.06  # percent/100
accounting_converterUnits.loc[
    idx["OMFix", "global", "CCGT", "2020"], "perUnitTotal"
] = 28.0  # Mio EUR per unit

accounting_converterUnits.loc[
    idx["Invest", "global", "OCGT", "2020"], "perUnitBuild"
] = 400.0  # Mio EUR per unit
accounting_converterUnits.loc[
    idx["Invest", "global", "OCGT", "2020"], "useAnnuity"
] = 1  # binary yes/no
accounting_converterUnits.loc[
    idx["Invest", "global", "OCGT", "2020"], "amorTime"
] = 30  # years
accounting_converterUnits.loc[
    idx["Invest", "global", "OCGT", "2020"], "interest"
] = 0.06  # percent/100
accounting_converterUnits.loc[
    idx["OMFix", "global", "OCGT", "2020"], "perUnitTotal"
] = 16.0  # Mio EUR per unit

accounting_converterUnits.loc[
    idx["Invest", "global", "LiIon", "2020"], "perUnitBuild"
] = 50.0  # Mio EUR per unit
accounting_converterUnits.loc[
    idx["Invest", "global", "LiIon", "2020"], "useAnnuity"
] = 1  # binary yes/no
accounting_converterUnits.loc[
    idx["Invest", "global", "LiIon", "2020"], "amorTime"
] = 25  # years
accounting_converterUnits.loc[
    idx["Invest", "global", "LiIon", "2020"], "interest"
] = 0.06  # percent/100
accounting_converterUnits.loc[
    idx["OMFix", "global", "LiIon", "2020"], "perUnitTotal"
] = 0.75  # Mio EUR per unit

accounting_converterUnits.loc[
    idx["Invest", "global", "PV", "2020"], "perUnitBuild"
] = 518.0
accounting_converterUnits.loc[idx["Invest", "global", "PV", "2020"], "useAnnuity"] = 1
accounting_converterUnits.loc[idx["Invest", "global", "PV", "2020"], "amorTime"] = 20
accounting_converterUnits.loc[idx["Invest", "global", "PV", "2020"], "interest"] = 0.06
accounting_converterUnits.loc[
    idx["OMFix", "global", "PV", "2020"], "perUnitTotal"
] = 7.7
accounting_converterUnits.loc[
    idx["Delta_PV", "global", "PV", "2020"], "perUnitBuild"
] = 1
accounting_converterUnits.loc[
    idx["Delta_PV", "global", "PV", "2020"], "perUnitDecom"
] = -1

accounting_converterUnits.loc[
    idx["Invest", "global", "GreenImport", "2040"], "perUnitBuild"
] = 500.0
accounting_converterUnits.loc[
    idx["Invest", "global", "GreenImport", "2040"], "useAnnuity"
] = 1
accounting_converterUnits.loc[
    idx["Invest", "global", "GreenImport", "2040"], "amorTime"
] = 30
accounting_converterUnits.loc[
    idx["Invest", "global", "GreenImport", "2040"], "interest"
] = 0.06
accounting_converterUnits.loc[
    idx["OMFix", "global", "GreenImport", "2040"], "perUnitTotal"
] = 25

accounting_converterUnits.loc[
    idx["Invest", "global", "WindOnshore", "2020"], "perUnitBuild"
] = 1173.0
accounting_converterUnits.loc[
    idx["Invest", "global", "WindOnshore", "2020"], "useAnnuity"
] = 1
accounting_converterUnits.loc[
    idx["Invest", "global", "WindOnshore", "2020"], "amorTime"
] = 25
accounting_converterUnits.loc[
    idx["Invest", "global", "WindOnshore", "2020"], "interest"
] = 0.06
accounting_converterUnits.loc[
    idx["OMFix", "global", "WindOnshore", "2020"], "perUnitTotal"
] = 25.8
accounting_converterUnits.loc[
    idx["Delta_Wind", "global", "WindOnshore", "2020"], "perUnitBuild"
] = 1
accounting_converterUnits.loc[
    idx["Delta_Wind", "global", "WindOnshore", "2020"], "perUnitDecom"
] = -1
accounting_converterUnits.dropna(how="all", inplace=True)
accounting_converterUnits.fillna(0, inplace=True)

m.parameter.add(accounting_converterUnits, "accounting_converterunits")
accounting_converterUnits
# %% [markdown]
# ### Sources and sinks
#
# #### Adding a demand profile as sink
#
# In this part, we set a demand for each node.
# The region names and years have to be included in the
# `map.aggregatenodesmodel` and `set.years` defined in the beginning.
# The name for the source-sink technology can be freely chosen.
#
# We need to specify that the demand is applied to the commodity "Elec" and
# that its profile needs to be matched exactly on an hour-by-hour level.

# %%
# "sourcesink_profile"
# load the profiles DataFrame, select the demand columns
demand = profiles[["demand_R1", "demand_R2", "demand_R3", "demand_R4"]]

demand = demand.div(1e3).mul(-1)
# transpose DataFrame for needed format
demand = demand.T

demand = demand.rename(
    index={
        "demand_R1": "R1_data",
        "demand_R2": "R2_data",
        "demand_R3": "R3_data",
        "demand_R4": "R4_data",
    }
)

demand_array = []
for year in m.set.years:
    current_demand = demand.copy()
    current_demand["years"] = year
    current_demand["techs"] = "Demand"
    current_demand["commodity"] = "Elec"
    current_demand["type"] = "fixed"
    current_demand = current_demand.set_index(
        ["years", "techs", "commodity", "type"], append=True
    )
    demand_array.append(current_demand)
all_demands = pd.concat(demand_array)
m.profile.add(all_demands, "sourcesink_profile")
all_demands.iloc[:, 0:8]
# %% [markdown]
# For the source-sink profile to be accurately understood by the model, we need
# to create a config with the information which profile and/or annual sum is
# actually used for which sources and sinks (in our case only "Demand").

# %%
# "sourcesink_config" (demand configuration)
sourcesink_config = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [m.set.nodesdata, m.set.yearssel, ["Demand"], ["Elec"]]
    )
)
sourcesink_config.loc[idx[m.set.nodesdata, :, :, :], "usesFixedProfile"] = 1
sourcesink_config.dropna(inplace=True)

m.parameter.add(sourcesink_config, "sourcesink_config")
sourcesink_config
# %% [markdown]
# #### Add sources for fuels and sinks for carbon emissions
#
# Since the gas turbine technologies use CH4 as a fuel, we need to allow import
# of CH4 for the model.
# This is very similar to the source-sink technology we used for the
# electrical demand.
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
sourcesink_annualSum.loc[idx[m.set.nodesdata, :, :, :], "upper"] = np.inf
sourcesink_annualSum.loc[idx[m.set.nodesdata, :, :, :], "lower"] = 0
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
sourcesink_config.loc[idx[m.set.nodesdata, :, :, :], "usesUpperSum"] = 1
sourcesink_config.loc[idx[m.set.nodesdata, :, :, :], "usesLowerProfile"] = 1
sourcesink_config.dropna(inplace=True)

m.parameter.add(sourcesink_config, "sourcesink_config")
sourcesink_config
# %%
# "accounting_sourcesinkFlow"
# setting a cost for methane imports
accounting_sourcesinkFlow = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["FuelCost", "Carbon"], ["global"], m.set.yearssel, ["FuelImport"], ["CH4"]]
    )
)
accounting_sourcesinkFlow["perFlow"] = 0.0
accounting_sourcesinkFlow.loc[
    idx["FuelCost", :, "2020", :], "perFlow"
] = 0.0306  # Mio EUR per GWh_ch CH4
accounting_sourcesinkFlow.loc[
    idx["FuelCost", :, "2030", :], "perFlow"
] = 0.0350  # Mio EUR per GWh_ch CH4
accounting_sourcesinkFlow.loc[
    idx["FuelCost", :, "2040", :], "perFlow"
] = 0.0400  # Mio EUR per GWh_ch CH4
accounting_sourcesinkFlow.loc[
    idx["FuelCost", :, "2050", :], "perFlow"
] = 0.0450  # Mio EUR per GWh_ch CH4

# Emission factor of Natural Gas: 0.2016 kton/GWh
# Calculated from: https://www.umweltbundesamt.de/sites/default/files/medien/479/publikationen/cc_28-2022_emissionsfaktoren-brennstoffe_bf.pdf
accounting_sourcesinkFlow.loc[idx["Carbon", :], "perFlow"] = 0.2016  # ktonnes/GWh

m.parameter.add(accounting_sourcesinkFlow, "accounting_sourcesinkflow")
accounting_sourcesinkFlow
# %% [markdown]
# ### Adding a transfer technology
#
# %% [markdown]
# We need to set up the transfer connections in the data by defining the
# starting and ending node of each link
# %%
# "transfer_linkStartEnd"
link_names = ["R1__R2", "R2__R3", "R1__R3", "R3__R4"]
data_nodes = m.set.nodesdata

transfer_connections = pd.DataFrame(
    index=pd.MultiIndex.from_product([link_names, data_nodes])
)
transfer_connections.loc[idx["R1__R2", "R1_data"], ["start"]] = 1
transfer_connections.loc[idx["R1__R2", "R2_data"], ["end"]] = 1
transfer_connections.loc[idx["R2__R3", "R2_data"], ["start"]] = 1
transfer_connections.loc[idx["R2__R3", "R3_data"], ["end"]] = 1
transfer_connections.loc[idx["R1__R3", "R1_data"], ["start"]] = 1
transfer_connections.loc[idx["R1__R3", "R3_data"], ["end"]] = 1
transfer_connections.loc[idx["R3__R4", "R3_data"], ["start"]] = 1
transfer_connections.loc[idx["R3__R4", "R4_data"], ["end"]] = 1
transfer_connections = transfer_connections.dropna(how="all").fillna(0)

m.parameter.add(transfer_connections, "transfer_linkstartend")
transfer_connections
# %% [markdown]
# Next we define the lengths for each corridor.
# We can use different link types to do so.
# %%
# "transfer_lengthParam"
link_types = ["land", "sea"]

transfer_lengths = pd.DataFrame(
    index=pd.MultiIndex.from_product([link_names, link_types])
)
transfer_lengths.loc[idx["R1__R2", "land"], ["length"]] = 1006.3
transfer_lengths.loc[idx["R2__R3", "land"], ["length"]] = 660.1
transfer_lengths.loc[idx["R1__R3", "land"], ["length"]] = 528.8
transfer_lengths.loc[idx["R3__R4", "land"], ["length"]] = 630.0
# transfer_lengths = transfer_lengths.dropna()

m.parameter.add(transfer_lengths, "transfer_lengthparam")
# %% [markdown]
# With the line corridors defined, we can start adding links to be optimized.
# %%
# "transfer_linksParam"
transfer_techs = ["HVDC"]

transfer_caps = read_dat("input/transfer_linksparam.dat")

m.parameter.add(transfer_caps, "transfer_linksparam")
transfer_caps
# %% [markdown]
# Define the technology information of the network
# %%
# "transfer_techParam"
tech_params = pd.DataFrame(index=pd.MultiIndex.from_product([transfer_techs, ["2020"]]))
tech_params.loc[:, "lifeTime"] = 40
tech_params.loc[:, "flowUpperLimit"] = 1

m.parameter.add(tech_params, "transfer_techparam")
tech_params
# %% [markdown]
# Define the commodity and rated capacity of the network technology
# %%
# "transfer_coefficient"
commodity = ["Elec"]

transfer_coefficient = pd.DataFrame(
    index=pd.MultiIndex.from_product([transfer_techs, m.set.yearssel, commodity])
)
transfer_coefficient["coefficient"] = 1  # GWh / h

m.parameter.add(transfer_coefficient, "transfer_coefficient")
transfer_coefficient
# %% [markdown]
# Define the losses for the converter stations
# %%
# "transfer_coefPerFlow"
coef_per_flow = pd.DataFrame(
    index=pd.MultiIndex.from_product([transfer_techs, ["2020"], commodity])
)
coef_per_flow[
    "coefPerFlow"
] = -0.014  # electrical losses of 14 MWh/h for each flow of 1 GWh/h

m.parameter.add(coef_per_flow, "transfer_coefperflow")
coef_per_flow

# %% [markdown]
# Define the losses for the links per km
# %%
# "transfer_coefPerLength"
coef_per_len = pd.DataFrame(
    index=pd.MultiIndex.from_product([transfer_techs, ["2020"], commodity, link_types])
)
coef_per_len.loc[
    idx[:, :, :, "land"], idx["coefPerLength"]
] = (
    -0.00004
)  # electrical losses of 40 kWh / h for each flow of 1 GWh / h and 1 km line length ~ 24 MWh / h for 600 km length
coef_per_len.loc[idx[:, :, :, "sea"], idx["coefPerLength"]] = -0.00003

m.parameter.add(coef_per_len, "transfer_coefperlength")
coef_per_len
# %% [markdown]
# Define indicators for each line built
# (for HVDC this is an AC/DC converter station at the beginning and end of the
# line)
# %%
# "accounting_transferLinks"
area = ["global"]

transfer_indicators = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["Invest", "OMFix"], area, transfer_techs, ["2020"]]
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
# %% [markdown]
# Define indicators for each line-km built (this needs the additional set for
# length-type modifiers, such as land and sea)
# %%
# "accounting_transferPerLength"
indicators_length = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["Invest", "OMFix", "OMVar"], area, transfer_techs, ["2020"], link_types]
    )
)
indicators_length.loc[idx["Invest", "global", :, :, "land"], "perLengthBuild"] = 0.2
indicators_length.loc[idx["Invest", "global", :, :, "land"], "interest"] = 0.06
indicators_length.loc[idx["Invest", "global", :, :, "land"], "amorTime"] = 40
indicators_length.loc[idx["Invest", "global", :, :, "land"], "useAnnuity"] = 1
indicators_length.loc[idx["OMFix", "global", :, :, "land"], "perLengthTotal"] = 0.004
indicators_length.loc[idx["OMVar", "global", :, :, "land"], "perFlow"] = 0.00002

indicators_length.loc[idx["Invest", "global", :, :, "sea"], "perLengthBuild"] = 0.4
indicators_length.loc[idx["Invest", "global", :, :, "sea"], "interest"] = 0.06
indicators_length.loc[idx["Invest", "global", :, :, "sea"], "amorTime"] = 40
indicators_length.loc[idx["Invest", "global", :, :, "sea"], "useAnnuity"] = 1
indicators_length.loc[idx["OMFix", "global", :, :, "sea"], "perLengthTotal"] = 0.005
indicators_length.loc[idx["OMVar", "global", :, :, "sea"], "perFlow"] = 0.00002
indicators_length = indicators_length.fillna(0)

m.parameter.add(indicators_length, "accounting_transferperlength")
indicators_length
# %% [markdown]
# #### The storage reservoir
#
# The storage features are always connected to a node and commodity combination
# and allow storing the connected commodity freely up to the rated capacity of
# the storage reservoir.
# We account for storage units in the same manner as for converter units and use
# a rated capacity to connect the units to a commodity and size.

# %%
# "storage_techParam"
storage_techParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([["LiIon"], ["2020"]])
)
storage_techParam.loc[idx["LiIon", :], "lifeTime"] = 25
storage_techParam.loc[idx["LiIon", :], "levelUpperLimit"] = 1

m.parameter.add(storage_techParam, "storage_techparam")
storage_techParam
# %% [markdown]
# For the storage size, we need to associate a commodity (here "Elec_LiIon") and
# a rated capacity for every storage reservoir unit.

# %%
# "storage_sizeParam"
# size of each storage unit
storage_sizeParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([["LiIon"], ["2020"], ["Elec_LiIon"]])
)
storage_sizeParam.loc[idx["LiIon", :, "Elec_LiIon"], "size"] = 8  # GWh_ch/unit
storage_sizeParam.loc[idx["LiIon", :, "Elec_LiIon"], "selfdischarge"] = -0.02
storage_sizeParam.dropna(inplace=True)

m.parameter.add(storage_sizeParam, "storage_sizeparam")
storage_sizeParam
# %% [markdown]
# We can set the storage reservoir upper limit to 30 units for a specific
# model region, therefore the model can build up to 240 GWh_ch of storage
# reservoir (8 GWh_ch / unit * 30 units = 240 GWh_ch).

# %%
# "storage_reservoirParam"
# installed storage reservoir units
storage_reservoirParam = read_dat("input/storage_reservoirparam.dat")

m.parameter.add(storage_reservoirParam, "storage_reservoirparam")
storage_reservoirParam
# %%
# "accounting_storageUnits"
# accounting for costs of storage
accounting_storageUnits = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["Invest", "OMFix"], ["global"], ["LiIon"], ["2020"]]
    )
)

accounting_storageUnits.loc[idx["Invest", :, :, :], "perUnitBuild"] = 200.0
accounting_storageUnits.loc[idx["Invest", :, :, :], "useAnnuity"] = 1
accounting_storageUnits.loc[idx["Invest", :, :, :], "amorTime"] = 25
accounting_storageUnits.loc[idx["Invest", :, :, :], "interest"] = 0.06
accounting_storageUnits.loc[idx["OMFix", :, :, :], "perUnitTotal"] = 3
accounting_storageUnits.fillna(0, inplace=True)

m.parameter.add(accounting_storageUnits, "accounting_storageunits")
accounting_storageUnits
# %% [markdown]
# ### Writing DataFrames to files
#
# %%
m.write(fileformat="dat")
# %% [markdown]
# We have finished building our REMix data model. In part b of this
# tutorial, we are looking into how we can execute it.
