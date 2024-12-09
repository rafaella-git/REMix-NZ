# %% [markdown]
# (tutorial_103_label)=
#
# # Tutorial 103 - Inter-regional energy transfer
#
# The image below shows the overview of different regions, that will be modeled
# in this tutorial.
# The individual regions are modeled identically to the previous tutrial, but
# they are connected with links.
#
# <div style="text-align: center;">
#
# ![Model regions overview for tutorial 103](../../img/REMix_tutorial103.svg "Model regions overview for tutorial 103")
#
# Model regions overview for tutorial 103
#
# </div>
#
# <div style="text-align: center;">
#
# ![Per region model for tutorial 103](../../img/REMix_tutorial102.svg "Per region model for tutorial 103")
#
# Per region model for tutorial 103
#
# </div>
#
# ## Part a: setting up the model
#
# In this tutorial we have a closer look at **transfer technologies**.
# In the second tutorial we had added technologies to store the electrical
# energy from the volatile renewable sources.
# As a next step, in this tutorial we include the possibility to transfer
# energy between the model nodes.
#
# As done before, we will use the previous tutorial_102 as a base model by
# reading its files into an Instance object `m` and adding a transfer
# technology to it.

# %%
# import dependencies
from remix.framework import Instance
import pandas as pd
import numpy as np
import pathlib as pt

# reading previous tutorial into Instance object `m`
_path_tut102_data = pt.Path("../tutorial_102/data")

if not _path_tut102_data.exists():
    raise IOError("You need to run tutorial 102a first!")

m = Instance.from_path(_path_tut102_data)

m.datadir = "./data"

# define often-used shortcut
idx = pd.IndexSlice
# %% [markdown]
# ### Adding a transfer technology
#
# After loading the model and dependencies from our base model, we can now
# simply add the components of the transfer technology.
# In this tutorial, this will be a high-voltage direct current ("HVDC") grid.

# %% [markdown]
# First we need to set up the transfer connections in the data by defining the
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
# We can use different link types.
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
transfer_lengths.loc[idx["R1__R2", "sea"], ["length"]] = 0.0
transfer_lengths.loc[idx["R2__R3", "sea"], ["length"]] = 0.0
transfer_lengths.loc[idx["R1__R3", "sea"], ["length"]] = 0.0
transfer_lengths.loc[idx["R3__R4", "sea"], ["length"]] = 0.0
# transfer_lengths = transfer_lengths.dropna()

m.parameter.add(transfer_lengths, "transfer_lengthparam")
# %% [markdown]
# With the  corridors now defined, we can start adding links to be optimized
# to the model.
# %%
# "transfer_linksParam"
transfer_techs = ["HVDC"]

transfer_caps = pd.DataFrame(
    index=pd.MultiIndex.from_product([link_names, m.set.yearssel, transfer_techs])
)
transfer_caps.loc[
    :, ["linksUpperLimit"]
] = 100  # Allow to build 100 GW for all links as the upper limit

m.parameter.add(transfer_caps, "transfer_linksparam")
transfer_caps
# %% [markdown]
# Define the technology information of the network
# %%
# "transfer_techParam"
tech_params = pd.DataFrame(
    index=pd.MultiIndex.from_product([transfer_techs, m.set.yearssel])
)
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
    index=pd.MultiIndex.from_product([transfer_techs, m.set.yearssel, commodity])
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
    index=pd.MultiIndex.from_product(
        [transfer_techs, m.set.yearssel, commodity, link_types]
    )
)
coef_per_len.loc[
    idx[:, :, :, "land"], idx["coefPerLength"]
] = (
    -0.00004
)  # electrical losses of 40 kWh / h for each flow of 1 GWh / h and 1 km  length ~ 24 MWh / h for 600 km length
coef_per_len.loc[idx[:, :, :, "sea"], idx["coefPerLength"]] = -0.00003

m.parameter.add(coef_per_len, "transfer_coefperlength")
coef_per_len
# %% [markdown]
# Define indicators for each  built
# (for HVDC this is an AC/DC converter station at the beginning and end of the
# )
# %%
# "accounting_transferLinks"
cost_indicators = ["Invest", "OMFix"]
area = ["global"]

transfer_indicators = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [cost_indicators, area, transfer_techs, m.set.yearssel]
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
# Define indicators for each -km built (this needs the additional set for
# length-type modifiers, such as land and sea)
# %%
# "accounting_transferPerLength"
indicators_length = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [cost_indicators, area, transfer_techs, m.set.yearssel, link_types]
    )
)
indicators_length.loc[idx["Invest", "global", :, :, "land"], "perLengthBuild"] = 0.544
indicators_length.loc[idx["Invest", "global", :, :, "land"], "interest"] = 0.06
indicators_length.loc[idx["Invest", "global", :, :, "land"], "amorTime"] = 40
indicators_length.loc[idx["Invest", "global", :, :, "land"], "useAnnuity"] = 1
indicators_length.loc[idx["OMFix", "global", :, :, "land"], "perLengthTotal"] = 0.00544

indicators_length.loc[idx["Invest", "global", :, :, "sea"], "perLengthBuild"] = 0.975
indicators_length.loc[idx["Invest", "global", :, :, "sea"], "interest"] = 0.06
indicators_length.loc[idx["Invest", "global", :, :, "sea"], "amorTime"] = 40
indicators_length.loc[idx["Invest", "global", :, :, "sea"], "useAnnuity"] = 1
indicators_length.loc[idx["OMFix", "global", :, :, "sea"], "perLengthTotal"] = 0.00975
indicators_length = indicators_length.fillna(0)

m.parameter.add(indicators_length, "accounting_transferperlength")
indicators_length
# %% [markdown]
# ### Adding additional demand and converters to the model
#
# Add profiles for PV and Wind onshore in for data nodes of R2, R3 and R4.
# %%
# "converter_activityProfile"
profiles = pd.read_csv("../_input/profiles.csv", index_col=0)

for data_node in ["R3_data", "R2_data", "R4_data"]:

    converter_activityProfile = profiles[["PV", "WindOnshore"]]

    # convert from MW to GW
    converter_activityProfile = converter_activityProfile.div(1e3).T

    converter_activityProfile = converter_activityProfile.div(
        converter_activityProfile.max(axis=1), axis=0
    )
    converter_activityProfile.index.names = ["techs"]

    converter_activityProfile["region"] = data_node
    converter_activityProfile["years"] = "2030"
    converter_activityProfile["type"] = "upper"
    converter_activityProfile = converter_activityProfile.reset_index().set_index(
        ["region", "years", "techs", "type"]
    )

    m.profile.add(converter_activityProfile, "converter_activityprofile")

m.profile.converter_activityprofile
converter_activityProfile.iloc[:, 0:8]
# %% [markdown]
# Add demand data nodes of R2, R3 and R4.
# %%
# "sourcesink_profile"
demand_R4_R2_CH = profiles[["demand_R4", "demand_R2", "demand_R3"]]

demand_R4_R2_CH = demand_R4_R2_CH.div(1e3).mul(-1)
# transpose DataFrame for needed format
demand_R4_R2_CH = demand_R4_R2_CH.T

demand_R4_R2_CH = demand_R4_R2_CH.rename(
    index={"demand_R4": "R4_data", "demand_R2": "R2_data", "demand_R3": "R3_data"}
)

# add columns and set them as index
demand_R4_R2_CH["years"] = "2030"
demand_R4_R2_CH["techs"] = "Demand"
demand_R4_R2_CH["commodity"] = "Elec"
demand_R4_R2_CH["type"] = "fixed"
demand_R4_R2_CH = demand_R4_R2_CH.set_index(
    ["years", "techs", "commodity", "type"], append=True
)

m.profile.add(demand_R4_R2_CH, "sourcesink_profile")
demand_R4_R2_CH.iloc[:, 0:8]
# %%
# "sourcesink_config" (demand configuration)
demand_cfg = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["R4_data", "R2_data", "R3_data"], m.set.yearssel, ["Demand"], ["Elec"]]
    )
)
demand_cfg["usesFixedProfile"] = 1

m.parameter.add(demand_cfg, "sourcesink_config")
demand_cfg
# %% [markdown]
# Add converter capacities to the data nodes of R2, R3 and R4.
# %%
# "converter_capacityParam"
converter_capacityParam = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [m.set.nodesdata, m.set.yearssel, ["CCGT", "PV", "WindOnshore"]]
    )
)
converter_capacityParam.loc[
    idx[("R2_data", "R3_data", "R4_data"), :, "CCGT"], "unitsUpperLimit"
] = 50  # GW_el
converter_capacityParam.loc[
    idx[("R2_data", "R3_data", "R4_data"), :, "CCGT"], "unitsLowerLimit"
] = 0  # GW_el
converter_capacityParam.loc[
    idx[("R2_data", "R3_data", "R4_data"), :, "WindOnshore"], "unitsUpperLimit"
] = 100  # GW_el
converter_capacityParam.loc[
    idx[("R2_data", "R3_data", "R4_data"), :, "WindOnshore"], "unitsLowerLimit"
] = 0  # GW_el
converter_capacityParam.loc[
    idx[("R2_data", "R3_data", "R4_data"), :, "PV"], "unitsUpperLimit"
] = 100  # GW_el
converter_capacityParam.loc[
    idx[("R2_data", "R3_data", "R4_data"), :, "PV"], "unitsLowerLimit"
] = 0  # GW_el

converter_capacityParam = converter_capacityParam.dropna(how="all")

m.parameter.add(converter_capacityParam, "converter_capacityparam")
converter_capacityParam
# %%
# limiting the annual sum of fuel imports into a model region
# "sourcesink_annualSum"
sourcesink_annualSum = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [m.set.nodesdata, m.set.yearssel, ["FuelImport"], ["CH4"]]
    )
)
sourcesink_annualSum.loc[
    idx[("R2_data", "R3_data", "R4_data"), :, :, :], "upper"
] = np.inf
sourcesink_annualSum.dropna(inplace=True)

m.parameter.add(sourcesink_annualSum, "sourcesink_annualsum")
sourcesink_annualSum

# %%
# limiting annual sum of carbon emissions
sourcesink_annualSum = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [m.set.nodesdata, m.set.yearssel, ["Emission"], ["CO2"]]
    )
)
sourcesink_annualSum.loc[
    idx[("R2_data", "R3_data", "R4_data"), :, :, :], "lower"
] = -np.inf
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
sourcesink_config.loc[
    idx[("R2_data", "R3_data", "R4_data"), :, :, :], "usesLowerSum"
] = 1
sourcesink_config.loc[
    idx[("R2_data", "R3_data", "R4_data"), :, :, :], "usesUpperProfile"
] = 1
sourcesink_config.dropna(inplace=True)

m.parameter.add(sourcesink_config, "sourcesink_config")
sourcesink_config

# %%
sourcesink_config = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [m.set.nodesdata, m.set.yearssel, ["FuelImport"], ["CH4"]]
    )
)
sourcesink_config.loc[
    idx[("R2_data", "R3_data", "R4_data"), :, :, :], "usesUpperSum"
] = 1
sourcesink_config.loc[
    idx[("R2_data", "R3_data", "R4_data"), :, :, :], "usesLowerProfile"
] = 1
sourcesink_config.dropna(inplace=True)

m.parameter.add(sourcesink_config, "sourcesink_config")
sourcesink_config
# %%
# writing files to `data/` directory
m.write(fileformat="dat")
# %% [markdown]
# That's it. We have successfully added tlinks to our model.
# We can now start a GAMS optimization run (part b).
