# %% [markdown]
# (tutorial_102_label)=
#
# # Tutorial 102 - Storage technologies
#
# <div style="text-align: center;">
#
# ![Model overview for tutorial 102](../../img/REMix_tutorial102.svg "Model overview for tutorial 102")
#
# Model overview of tutorial 102
#
# </div>
#
# ## Part a: setting up the model
#
# In this tutorial we have a closer look at **storage technologies**.
# In the first tutorial we had renewable energies in the system and checked two
# weeks with the highest and lowest renewable generation.
# The feed-in from renewable energies was mainly limited by the feed-in
# profiles.
# As a next step, in this tutorial we include technologies to store the
# electrical energy from the volatile renewable sources and thus add a
# flexibility.
#
# As mentioned during tutorial_101, we will use it as a base model here by
# reading its files into an Instance object `m` and adding a storage technology
# to it.

# %%
# importing dependencies
from remix.framework.api.instance import Instance
import pandas as pd
import pathlib as pt

# reading in model built in `tutorial_101a_build`
_path_tut1_data = pt.Path("../tutorial_101/data")

if not _path_tut1_data.exists():
    raise IOError("You need to run tutorial 1a first!")

m = Instance.from_path(_path_tut1_data)

m.datadir = "./data"

# define often-used shortcut
idx = pd.IndexSlice
# %% [markdown]
# ### Adding a storage technology
#
# After loading the model and dependencies from our base model (i.e.
# `tutorial_101a_build.py`), we can now simply add the components of the
# storage.
#
# Storage technologies are typically comprised of two parts:
# (1) the energy storage itself;
# (2) the component for charging and discharging the storage.
#
# Similarly, in REMix the storages are also built on top of two different
# components.
# A storage converter for charging and discharging a storage reservoir and that
# reservoir itself that contains the chosen commodity (in this case
# electricity).
#
# #### The charging/discharging unit (=converter)
#
# First, we will define the storage converter, i.e. the charging/discharging
# unit.
#
# We can use the same features we used for the converters of conventional power
# plants.
# The difference is that a storage by definition converts one commodity into the
# same commodity (e.g. electricity to electricity).
#
# As an example, we introduce a lithium-ion battery as electricity storage.

# %%
# "converter_techParam"
converter_techParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([["Battery"], m.set.yearssel])
)
converter_techParam.loc[idx["Battery", :], "lifeTime"] = 25
converter_techParam.loc[idx["Battery", :], "activityUpperLimit"] = 1

m.parameter.add(converter_techParam, "converter_techparam")
converter_techParam
# %%
# "converter_capacityParam"
converter_capacityParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([m.set.nodesdata, m.set.yearssel, ["Battery"]])
)
converter_capacityParam.loc[
    idx["R1_data", :, "Battery"], "unitsUpperLimit"
] = 30  # GW_el
converter_capacityParam.dropna(inplace=True)

m.parameter.add(converter_capacityParam, "converter_capacityparam")
converter_capacityParam
# %% [markdown]
# In contrast to the previous modeling of converter units for conventional power plants, we now need to define a
# reversible activity. In this example, we can both charge and discharge our lithium-ion battery with the same power
# unit. Therefore, we add both activities---`Charge` and `Discharge`---and use the coefficients to model the
# corresponding losses.
#
# We can also use two different converters for charging and discharging. This is necessary when wanting to better
# represent the real-world difference between the turbine and optional pumps in hydroelectric power plants for example.
# These can then also have different rated powers.
#
# A storage in REMix per definition has the same input and output commodity. To be able to account for storage losses,
# it is necessary to define a dummy commodity (here called `Elec_LiIon`), which is only used inside that one technology.
#
# In this tutorial, we fill the two activities of our single converter unit for charging and discharging so that each
# process has an efficiency of 95 %.

# %%
# "converter_coefficient"
converter_coefficient = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["Battery"], m.set.yearssel, ["Charge", "Discharge"], ["Elec", "Elec_LiIon"]]
    )
)

converter_coefficient.loc[
    idx["Battery", :, "Charge", "Elec"], "coefficient"
] = -1  # GW_el
converter_coefficient.loc[
    idx["Battery", :, "Charge", "Elec_LiIon"], "coefficient"
] = 0.95  # GW_el in LiIon
converter_coefficient.loc[
    idx["Battery", :, "Discharge", "Elec"], "coefficient"
] = 1  # GW_el
converter_coefficient.loc[
    idx["Battery", :, "Discharge", "Elec_LiIon"], "coefficient"
] = -1.05  # GW_el in LiIon
converter_coefficient.dropna(how="all", inplace=True)

m.parameter.add(converter_coefficient, "converter_coefficient")
converter_coefficient
# %%
# "accounting_converterUnits"
accounting_converterUnits = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["Invest", "OMFix"], ["global"], ["Battery"], m.set.yearssel]
    )
)

accounting_converterUnits.loc[
    idx["Invest", "global", "Battery", "2030"], "perUnitBuild"
] = 50  # million EUR / unit
accounting_converterUnits.loc[
    idx["Invest", "global", "Battery", "2030"], "useAnnuity"
] = 1  # binary yes/no
accounting_converterUnits.loc[
    idx["Invest", "global", "Battery", "2030"], "amorTime"
] = 25  # years
accounting_converterUnits.loc[
    idx["Invest", "global", "Battery", "2030"], "interest"
] = 0.06  # percent/100
accounting_converterUnits.loc[
    idx["OMFix", "global", "Battery", "2030"], "perUnitTotal"
] = 0.75  # million EUR per unit and year
accounting_converterUnits.fillna(0, inplace=True)

m.parameter.add(accounting_converterUnits, "accounting_converterunits")
accounting_converterUnits
# %% [markdown]
# #### The storage reservoir
#
# The storage features are always connected to a node and commodity combination
# and allow storing the connected commodity freely up to the rated capacity of
# the storage reservoir.
# We account for storage units in the same manner as for converter units and use
# a rated capacity to connect the units to a commodity and size.
# Storage technologies and converter technologies have the same name to make it
# easier to represent them as the same technology.

# %%
# "storage_techParam"
storage_techParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([["Battery"], m.set.yearssel])
)
storage_techParam.loc[idx["Battery", :], "lifeTime"] = 25
storage_techParam.loc[idx["Battery", :], "levelUpperLimit"] = 1

m.parameter.add(storage_techParam, "storage_techparam")
storage_techParam
# %% [markdown]
# For the storage size, we need to associate a commodity (here "Elec_LiIon") and
# a rated capacity for every storage reservoir unit.

# %%
# "storage_sizeParam"
# size of each storage unit
storage_sizeParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([["Battery"], m.set.yearssel, ["Elec_LiIon"]])
)
storage_sizeParam.loc[idx["Battery", :, "Elec_LiIon"], "size"] = 8  # GWh_ch/unit
storage_sizeParam.dropna(inplace=True)

m.parameter.add(storage_sizeParam, "storage_sizeparam")
storage_sizeParam
# %% [markdown]
# Now we can set the storage reservoir upper limit to 30 units for a specific
# model region, therefore the model can build up to 240 GWh_ch of storage
# reservoir (8 GWh_ch / unit * 30 units = 240 GWh_ch).

# %%
# "storage_reservoirParam"
# installed storage reservoir units
storage_reservoirParam = pd.DataFrame(
    index=pd.MultiIndex.from_product([m.set.nodesdata, m.set.yearssel, ["Battery"]])
)
storage_reservoirParam.loc[
    idx["R1_data", :, "Battery"], "unitsUpperLimit"
] = 30  # units
storage_reservoirParam.dropna(inplace=True)

m.parameter.add(storage_reservoirParam, "storage_reservoirparam")
storage_reservoirParam
# %%
# "accounting_storageUnits"
# accounting for costs of storage
accounting_storageUnits = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [["Invest", "OMFix"], ["global"], ["Battery"], m.set.yearssel]
    )
)

accounting_storageUnits.loc[idx["Invest", :, :, :], "perUnitBuild"] = (
    105.5 * 8
)  # Since our storage unit can store 8 GWh we need to scale the million EUR/GWh value with 8
accounting_storageUnits.loc[idx["Invest", :, :, :], "useAnnuity"] = 1
accounting_storageUnits.loc[idx["Invest", :, :, :], "amorTime"] = 25
accounting_storageUnits.loc[idx["Invest", :, :, :], "interest"] = 0.06
accounting_storageUnits.loc[idx["OMFix", :, :, :], "perUnitTotal"] = 105.5 * 8 * 0.015
accounting_storageUnits.fillna(0, inplace=True)

m.parameter.add(accounting_storageUnits, "accounting_storageunits")
accounting_storageUnits
# %%
# write all files to `data/` directory
m.write(fileformat="dat")
# %% [markdown]
# That's it. We have successfully added a lithium-ion battery as storage
# technology to our model. We can now start a GAMS optimization run (part b).
