# %% [markdown]
# ## Part c: evaluation of results
#
# %%
# importing dependencies
from remix.framework.tools.gdx import GDXEval
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

technology_colors = {
    "CCGT": "#ff5f2d",
    "OCGT": "#e7864b",
    "PV": "#ffc000",
    "WindOnshore": "#9dc3e6",
    "HVDC": "#70ad47",
    "LiIon": "#cc99ff",
}
background_color = "#fffaebff"

result_dir = "./results"

# define often-used shortcut
idx = pd.IndexSlice
# %%
# read in the output `*.gdx` file from the optimization in GAMS
results_loose = GDXEval(f"{result_dir}/tutorial_202_loose.gdx")
results_tight = GDXEval(f"{result_dir}/tutorial_202_tight.gdx")
results_infinite = GDXEval(f"{result_dir}/tutorial_202_infinite.gdx")
# %% [markdown]
# ### Evaluating converter capacities

# %%
# convert converter capacities to a Pandas DataFrame
caps = results_loose["converter_caps"]
caps = caps[caps > 0.01].dropna()  # Remove all capacities with less than 10 MW

caps.loc[idx[:, "2050", :, "Elec", "total"], :].round(2)
# %% [markdown]
# We can now check the installed connection capacities between the different
# model nodes. This is done in a similar fashion as with the converter
# capacities.
# %%
transfer_caps = results_loose["transfer_caps"]

transfer_caps.loc[idx[:, :, :, "2050", :, "Elec", "total"], :].round(2)
# %% [markdown]
#
# To get a feeling on where we benefit from the electrical network, we can check
# the annual transferred energy between model nodes. Since the network contains
# information on the direction of flows, we need to also account for the
# direction. The energy flow from model region A (nodesModel) to B
# (nodesModel_a) is defined as positive, whereas the flow from B to A is
# accounted negative. With the "balanceType" entry we can check for the
# individual flows from A to B (positive), flows from B to A (negative), annual
# sum of directed flows (netto = positive + negative), annual sum of energy
# transferred (brutto = positive - negative), or line utilization
# (flh = brutto / line capacity).

# %%
transfer_flows = results_loose["transfer_flows_annual"]
transfer_flows = transfer_flows[
    transfer_flows.abs() > 0.1
].dropna()  # Remove all flows with less than 0.1 GWh

transfer_flows = (
    transfer_flows.loc[idx[:, :, :, :, :, "Elec", "netto"], :].div(1e3).round(2)
)  # Convert to TWh

transfer_flows
# %% [markdown]
#  We identified `R1_model` as the main importing model region, so we can
# further check the behavior of the hourly electricity supply.

# %%
# visualization of commodity balance for model region R1_model
commodities = results_loose["commodity_balance"]


def plot_highest_transmission(commodities, name="Loose", region="R1_model"):
    elec_R1 = (
        commodities.loc[idx[:, region, "2050", :, "Elec"], :]
        .groupby(["timeModel", "techs"])
        .sum()
        .unstack("techs")
    )

    demand_R1 = elec_R1.loc[:, idx[:, "Demand"]]
    demand_R1.columns = demand_R1.columns.get_level_values(1)

    positive_R1 = elec_R1.drop(columns=("value", "Demand"))
    positive_R1 = positive_R1[positive_R1 > 0].fillna(0)
    positive_R1.columns = positive_R1.columns.get_level_values(1)

    negative_R1 = elec_R1.drop(columns=("value", "Demand"))
    negative_R1 = negative_R1[negative_R1 < 0].fillna(0)
    negative_R1.columns = negative_R1.columns.get_level_values(1)

    hours_per_interval = 168  # 168 hours per week
    rolling_mean = positive_R1[["HVDC"]].sum(axis=1).rolling(hours_per_interval).mean()
    last_hour_max = rolling_mean.argmax()

    timeslice = range(last_hour_max - hours_per_interval, last_hour_max)

    fig, ax1 = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor(background_color)
    ax1.set_facecolor(background_color)
    positive_R1.iloc[timeslice].plot.area(stacked=True, ax=ax1, color=technology_colors)

    plt.legend(loc="upper left")
    plt.ylabel("Generation in GWh_el")
    plt.title(
        f"Generation in week with highest transmission to model region R1: {name}"
    )
    plt.ylim(bottom=0, top=150)
    ax2 = ax1.twinx()
    demand_R1.iloc[timeslice].mul(-1).plot(
        kind="line", stacked=True, ax=ax2, color="black"
    )
    plt.legend(loc="upper right")
    plt.ylabel("Demand in GWh_el")
    ax2.set_ylim(ax1.get_ylim())
    plt.ylim(bottom=0, top=150)

    fig.subplots_adjust(bottom=0.1 * demand_R1.index.nlevels)


plot_highest_transmission(results_infinite["commodity_balance"], "Infinite")
plot_highest_transmission(results_loose["commodity_balance"], "Loose")
plot_highest_transmission(results_tight["commodity_balance"], "Tight")
# %%
# map defaults and options
shp_file = "../_input/shp_file/tutorial"
shp_attrcol = "ISO"

lat = [35.0, 57.0]
lon = [-6.0, 20.0]

centroids = {
    "FR": (2.34, 47.15),
    "DE": (10.18, 51.35),
    "IT": (12.80, 42.61),
    "CH": (8.17, 46.86),
}

plt.rcParams.update({"figure.autolayout": True})  # use tight_layout
plt.rcParams.update({"figure.titlesize": 20})  # size subtitle
plt.rcParams.update({"figure.dpi": 75})  # size subtitle
plt.rcParams.update({"savefig.dpi": 300})  # dpi
plt.rcParams.update({"font.size": 16})  # default font size
plt.rcParams.update({"axes.titlesize": 20})  # title size per subplot
plt.rcParams.update({"axes.labelsize": 18})  # label size colormap

figsize_dual = (13.0, 6)
figsize_single = (7.5, 6)

# %%
# visualization of energy generation and flows between model regions
def plot_generation_flows(df_annual):
    df_annual = (
        df_annual.reset_index()
        # .replace({"R3_model": "CH", "R1_model": "DE", "R2_model": "FR", "R4_model": "IT"})
        .set_index(df_annual.index.names)
    )
    df_annual = df_annual.loc[
        idx[
            ["R1_model", "R2_model", "R3_model", "R4_model"], "2050", :, "Elec", "netto"
        ]
    ]

    map_pv = (
        df_annual.loc[idx[:, :, "PV"], idx[:]].groupby("accNodesModel").sum().div(1e3)
    )
    map_wind = (
        df_annual.loc[idx[:, :, "WindOnshore"], idx[:]]
        .groupby("accNodesModel")
        .sum()
        .div(1e3)
    )

    g = nx.Graph()
    g.add_nodes_from(
        [
            (
                node,
                {"pv": map_pv.loc[node, "value"], "wind": map_wind.loc[node, "value"]},
            )
            for node in map_wind.index.values
        ]
    )
    n2n_flow = results_loose["transfer_flows_annual"]
    n2n_flow = (
        n2n_flow.reset_index()
        # .replace({"R3_model": "CH", "R1_model": "DE", "R2_model": "FR", "R4_model": "IT"})
        .set_index(n2n_flow.index.names)
    )

    n2n_flow = n2n_flow[n2n_flow != 0].dropna(how="all")

    flow_elec = (
        n2n_flow[n2n_flow != 0]
        .loc[idx[:, :, :, :, :, "Elec", "netto"]]
        .dropna()
        .groupby(["nodesModel_start", "nodesModel_end"])
        .sum()
        .div(1e3)
        .abs()
    )

    flow_elec = flow_elec[flow_elec > 0.1].rename(columns={"value": "flow_elec"})
    flow_flh = (
        n2n_flow[n2n_flow != 0]
        .loc[idx[:, :, :, :, :, "Elec", "flh"]]
        .dropna()
        .groupby(["nodesModel_start", "nodesModel_end"])
        .sum()
        .abs()
    )

    flow_flh = flow_flh[flow_flh > 10].rename(columns={"value": "flh"})

    flow_edges = pd.concat([flow_flh, flow_elec], axis=1).dropna().reset_index()
    for _, row in flow_edges.iterrows():
        g.add_edge(
            row["nodesModel_start"],
            row["nodesModel_end"],
            flow_elec=row["flow_elec"],
            flh=row["flh"],
        )
    pos = nx.circular_layout(g)
    fig = plt.figure(figsize=figsize_dual)
    fig.patch.set_facecolor(background_color)
    ax1 = fig.add_subplot(121)

    nx.draw(
        g,
        pos,
        ax=ax1,
        with_labels=True,
        node_color=technology_colors["WindOnshore"],
        width=[0.2 * g[u][v]["flow_elec"] for u, v in g.edges()],
        node_size=[25 * g.nodes[n]["wind"] for n in g.nodes],
        edge_color=technology_colors["HVDC"],
    )
    ax1.set_facecolor(background_color)
    plt.title("WindOnshore")
    ax1.axis("off")
    ax2 = fig.add_subplot(122)
    nx.draw(
        g,
        pos,
        ax=ax2,
        node_color=technology_colors["PV"],
        with_labels=True,
        width=[0.2 * g[u][v]["flow_elec"] for u, v in g.edges()],
        node_size=[25 * g.nodes[n]["pv"] for n in g.nodes],
        edge_color=technology_colors["HVDC"],
    )
    ax2.set_facecolor(background_color)
    plt.title("PV")
    ax2.axis("off")
    fig.set_facecolor(background_color)

    return fig


# Infinite Carbon Budget
fig1 = plot_generation_flows(results_infinite["commodity_balance_annual"])
fig1.suptitle("Budget: Infinite")
# Loose Carbon Budget
fig2 = plot_generation_flows(results_loose["commodity_balance_annual"])
fig2.suptitle("Budget: Loose")
# Strict Carbon Budget
fig3 = plot_generation_flows(results_tight["commodity_balance_annual"])
fig3.suptitle("Budget: Tight")
# %%
# total capacity per year / annual energy per year / carbon emissions per year
gen_capacities = results_loose["converter_caps"]

years = ["2020", "2030", "2040", "2050"]


def plot_total_capacities(df, ax, nodesData, techs, years, upper=800):
    nodesData = [n for n in nodesData if n in df.index.get_level_values(0).unique()]
    years = [n for n in years if n in df.index.get_level_values(1).unique()]
    techs = [n for n in techs if n in df.index.get_level_values(2).unique()]
    generation_capacities = df.loc[idx[nodesData, years, techs, :, "total"]]
    generation_capacities = (
        generation_capacities.reset_index()
        .groupby(["accNodesModel", "accYears", "techs"])
        .agg({"value": sum})
        .loc[idx[nodesData, years, techs], :]
    )
    if len(techs) == 1:
        generation_capacities = generation_capacities.droplevel("techs")
        generation_capacities = generation_capacities.unstack()
        generation_capacities.columns = generation_capacities.columns.get_level_values(
            1
        )
        generation_capacities.plot.bar(
            stacked=False, ax=ax, color=["#cd7f32", "#b7410e", "#465945", "#7ba05b"]
        )
        plt.xticks(rotation=0)
    else:
        generation_capacities = generation_capacities.unstack()
        generation_capacities.columns = generation_capacities.columns.get_level_values(
            1
        )
        generation_capacities.plot.bar(stacked=True, ax=ax, color=technology_colors)
        # plt.xticks(rotation=0)
    plt.ylim(bottom=0, top=upper)


# generation capacities
generation_techs = ["CCGT", "OCGT", "PV", "WindOnshore"]
nodes = ["R1_model", "R2_model", "R3_model", "R4_model"]
# Infinite
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
plot_total_capacities(
    results_infinite["converter_caps"], ax1, nodes, generation_techs, years
)
plt.legend(loc="upper left")
plt.ylabel("Total capacity in GW")
plt.title("Installed generator capacities of R1_model: Infinite")
# Loose
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
plot_total_capacities(
    results_loose["converter_caps"], ax1, nodes, generation_techs, years
)
plt.legend(loc="upper left")
plt.ylabel("Total capacity in GW")
plt.title("Installed generator capacities of R1_model: Loose")
# Tight
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
plot_total_capacities(
    results_tight["converter_caps"], ax1, nodes, generation_techs, years
)
plt.legend(loc="upper left")
plt.ylabel("Total capacity in GW")
plt.title("Installed generator capacities: Tight")
# %%
# storage capacities
storage_techs = ["LiIon"]
# Infinite
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
plot_total_capacities(
    results_infinite["storage_caps"], ax1, nodes, storage_techs, years, upper=1500
)
plt.legend(loc="upper left")
plt.ylabel("Total storage capacity in GWh")
plt.title("Installed storage capacities: Infinite")
# Loose
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
plot_total_capacities(
    results_loose["storage_caps"], ax1, nodes, storage_techs, years, upper=1500
)
plt.legend(loc="upper left")
plt.ylabel("Total storage capacity in GWh")
plt.title("Installed storage capacities: Loose")
# Tight
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
plot_total_capacities(
    results_tight["storage_caps"], ax1, nodes, storage_techs, years, upper=1500
)
plt.legend(loc="upper left")
plt.ylabel("Total storage capacity in GWh")
plt.title("Installed storage capacities: Tight")
years = ["2020", "2030", "2040", "2050"]

# %%
def plot_year_generation(df, ax, nodesData, techs, commodity, years, upper=600):
    generation_capacities = df.loc[idx[nodesData, years, techs, commodity, "netto"]]
    generation_capacities = (
        generation_capacities.reset_index()
        .astype({col: object for col in generation_capacities.index.names})
        .groupby(["accNodesModel", "accYears", "techs"])
        .agg({"value": sum})
    )
    if len(generation_capacities.index.levels[0]) == 1:
        generation_capacities = generation_capacities.groupby(
            ["accYears", "techs"]
        ).agg({"value": sum})
    generation_capacities = generation_capacities.unstack() / 1000
    generation_capacities.columns = generation_capacities.columns.get_level_values(1)

    generation_capacities.plot.bar(stacked=True, ax=ax, color=technology_colors)
    # plt.xticks(rotation=0)
    plt.ylim(bottom=0, top=upper)


# generation capacities
generation_techs = ["CCGT", "OCGT", "PV", "WindOnshore"]
# Infinite
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
plot_year_generation(
    results_infinite["commodity_balance_annual"],
    ax1,
    nodes,
    generation_techs,
    "Elec",
    years,
)
plt.legend(bbox_to_anchor=(1.0, 0.70))
plt.ylabel("Total annual generation in TWh_el")
plt.title("Annual generation: Infinite")
# Loose
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
plot_year_generation(
    results_loose["commodity_balance_annual"],
    ax1,
    nodes,
    generation_techs,
    "Elec",
    years,
)
plt.legend(bbox_to_anchor=(1.0, 0.70))
plt.ylabel("Total annual generation in TWh_el")
plt.title("Annual generation: Loose")
# Tight
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
plot_year_generation(
    results_tight["commodity_balance_annual"],
    ax1,
    nodes,
    generation_techs,
    "Elec",
    years,
)
plt.legend(bbox_to_anchor=(1.0, 0.70))
plt.ylabel("Total annual generation in TWh_el")
plt.title("Annual generation: Tight")
# %%
# Infinite
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
plot_year_generation(
    results_infinite["commodity_balance_annual"],
    ax1,
    "global",
    generation_techs,
    "Elec",
    years,
    upper=1500,
)
plt.legend(bbox_to_anchor=(1.0, 0.70))
plt.ylabel("Total annual generation in TWh_el")
plt.title("Annual generation of the whole model: Infinite")
# Loose
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
plot_year_generation(
    results_loose["commodity_balance_annual"],
    ax1,
    "global",
    generation_techs,
    "Elec",
    years,
    upper=1500,
)
plt.legend(bbox_to_anchor=(1.0, 0.70))
plt.ylabel("Total annual generation in TWh_el")
plt.title("Annual generation of the whole model: Loose")
# Tight
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
plot_year_generation(
    results_tight["commodity_balance_annual"],
    ax1,
    "global",
    generation_techs,
    "Elec",
    years,
    upper=1500,
)
plt.legend(bbox_to_anchor=(1.0, 0.70))
plt.ylabel("Total annual generation in TWh_el")
plt.title("Annual generation of the whole model: Tight")
# %%
indicator_accounting = results_loose["indicator_accounting_comp"]


def plot_emissions_by_node(df, ax, nodesModel, years, indicator):
    emissions = df.loc[
        idx["global", "horizon", indicator, nodesModel, years, indicator], "value"
    ]
    emissions = (
        emissions.reset_index()
        .astype({col: object for col in emissions.index.names})
        .groupby(["accNodesModel_a", "accYears_a"])
        .agg({"value": sum})
        / 1000
    )["value"].unstack()
    emissions.index.name = "Accounting Node"
    emissions.columns.name = "Accounting Year"
    emissions.plot.bar(ax=ax, color=["#cd7f32", "#b7410e", "#465945", "#7ba05b"])
    plt.ylim(bottom=0, top=1800)


# Infinite
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
plot_emissions_by_node(
    results_infinite["indicator_accounting_comp"],
    ax1,
    ["R1_model", "R2_model", "R3_model", "R4_model"],
    ["2020", "2030", "2040", "2050"],
    "Carbon",
)
plt.xticks(rotation=0)
plt.legend(bbox_to_anchor=(1.0, 0.70))
plt.ylabel("CO2 equivalent in Mtonnes")
plt.title("Annual emissions: Infinite")
# Loose
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
plot_emissions_by_node(
    results_loose["indicator_accounting_comp"],
    ax1,
    ["R1_model", "R2_model", "R3_model", "R4_model"],
    ["2020", "2030", "2040", "2050"],
    "Carbon",
)
plt.xticks(rotation=0)
plt.legend(bbox_to_anchor=(1.0, 0.70))
plt.ylabel("CO2 equivalent in Mtonnes")
plt.title("Annual emissions: Loose")
# Tight
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
plot_emissions_by_node(
    results_tight["indicator_accounting_comp"],
    ax1,
    ["R1_model", "R2_model", "R3_model", "R4_model"],
    ["2020", "2030", "2040", "2050"],
    "Carbon",
)
plt.xticks(rotation=0)
plt.legend(bbox_to_anchor=(1.0, 0.70))
plt.ylabel("CO2 equivalent in Mtonnes")
plt.title("Annual emissions: Tight")

# %%
indicator = results_loose["indicator_accounting"]
indicator_tight = results_tight["indicator_accounting"]
indicator_infinite = results_infinite["indicator_accounting"]

fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)


def plot_emissions_by_budget(df0, df1, df2, ax):
    df = (
        (
            pd.concat(
                {
                    "tight": filter_emissions(df0),
                    "loose": filter_emissions(df1),
                    "infinite": filter_emissions(df2),
                }
            )
            / 100
        )
        .swaplevel(0, 1)
        .unstack()
    )
    df.plot(
        ax=ax, color=["#70ad47", "#cc99ff", "#b7410e"]
    )  # , linestyle=["-.", "-", "."])
    plt.ylabel("CO2 equivalent in Mtonnes")
    plt.title("Annual system emissions")
    return df


def filter_emissions(df):
    df = (
        df.loc[
            idx[
                ["R1_model", "R2_model", "R3_model", "R4_model"],
                ["2020", "2030", "2040", "2050"],
                "Carbon",
            ],
            "value",
        ]
        .groupby("accYears")
        .sum()
    )
    return df


df = plot_emissions_by_budget(indicator_tight, indicator, indicator_infinite, ax1)
plt.show()
