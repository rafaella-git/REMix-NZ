# %%
import os
import time
import yaml
import fiona
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path

# from shapely.geos import lgeos
from shapely.geometry import Polygon, MultiPolygon, LineString, base, shape
from mpl_toolkits.basemap import Basemap
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon as PolyPatch
from matplotlib.collections import PatchCollection, LineCollection
from remix.framework.tools.gdx import GDXEval
import plotly.express as px

sns.set(color_codes=True)
sns.set_style("whitegrid")
idx = pd.IndexSlice

plt.rcParams.update({"figure.titlesize": 20})  # size suptitle
plt.rcParams.update({"figure.dpi": 75})  # size suptitle
plt.rcParams.update({"savefig.dpi": 300})  # dpi
plt.rcParams.update({"font.size": 16})  # default font size
plt.rcParams.update({"axes.titlesize": 20})  # title size per subplot
plt.rcParams.update({"axes.labelsize": 18})  # label size colormap
plt.rcParams.update({"xtick.labelsize": 16})
plt.rcParams.update({"ytick.labelsize": 16})
plt.rcParams.update({"legend.fontsize": 14})
# plt.rcParams.keys()

figsize_dual = (13.4, 7)
figsize_single = (9, 7)

# %%%

def make_valid(ob):
    return ob if ob.is_valid else base.geom_factory(lgeos.GEOSMakeValid(ob._geom))

# %% 
def redistribute_vertices(geom, distance):
    if geom.geom_type == 'LineString':
        num_vert = int(round(geom.length / distance))
        if num_vert == 0:
            num_vert = 1
        return LineString(
            [geom.interpolate(float(n) / num_vert, normalized=True)
             for n in range(num_vert + 1)])
    elif geom.geom_type == 'MultiLineString':
        parts = [redistribute_vertices(part, distance)
                 for part in geom]
        return type(geom)([p for p in parts if not p.is_empty])
    else:
        raise ValueError('unhandled geometry %s', (geom.geom_type,))


def load_polypatches(shp, shp_attr):
    coll = {}
    with fiona.open(shp, "r") as f:
        for i in f:
            coll[i["properties"][shp_attr]] = make_valid(shape(i["geometry"]))

    patch_coll = {}
    for k, i in coll.items():
        if isinstance(i, MultiPolygon):
            patch_coll[k] = [PolyPatch(np.array(j.exterior.coords.xy).T, True) for j in i.geoms]
        elif isinstance(i, Polygon):
            patch_coll[k] = [PolyPatch(np.array(i.exterior.coords.xy).T, True)]
        elif isinstance(i, LineString):
            patch_coll[k] = i # [np.array(i.coords.xy).T]
            if "__" in k:
                k_r = "__".join(reversed(k.split("__")))
                patch_coll[k_r] = i.reverse() # [np.fliplr(np.array(i.coords.xy)).T]

    return patch_coll


def add_basemap(lat, lon):
    m = Basemap(
        llcrnrlon=lon[0],
        urcrnrlon=lon[1],
        llcrnrlat=lat[0],
        urcrnrlat=lat[1],
        resolution="i",
        projection="cyl",
        fix_aspect=True,
    )

    m.fillcontinents(color="#e0e0e0")


def get_cmap(cmap, data, limits):
    cmap = matplotlib.cm.get_cmap(cmap)
    if limits is None:
        limits = [None, None]
    if limits[0] is None:
        limits[0] = float(data.min())
    if limits[1] is None:
        limits[1] = float(data.max())
    return cmap, limits


def plot_lines(ax, data, patch_coll, cmap, limits):
    for ikey in data.index.values:
        if ikey in patch_coll:
            ax.add_collection(
                LineCollection(
                    [np.array(patch_coll[ikey].coords.xy).T],
                    edgecolor=cmap(float((data.loc[ikey] - limits[0]) / (limits[1] - limits[0]))),
                    linewidths=8 * float((data.loc[ikey] - limits[0]) / (limits[1] - limits[0])),
                    zorder=2,
                )
            )
        else:
            print(f"Dataframe index {ikey} has no corresponding shape")


def plot_arrows(ax, data, patch_coll, cmap, limits):
    for ikey in data.index.values:
        if ikey in patch_coll and ikey in data.index:
                ls = LineString(patch_coll[ikey])

                mls = redistribute_vertices(ls, 1)
                x_position = np.array([i for i in mls.coords.xy[0]])
                y_position = np.array([i for i in mls.coords.xy[1]])
                x_direction = np.zeros(shape=x_position.shape).round(3)
                y_direction = np.zeros(shape=y_position.shape).round(3)
                x_direction[1:-1] = (x_position[2:] - x_position[:-2])
                y_direction[1:-1] = (y_position[2:] - y_position[:-2])
                arrow_len = np.sqrt(np.square(x_direction) + np.square(y_direction)) + .001
                x_direction = x_direction / arrow_len
                y_direction = y_direction / arrow_len

                scale = 5 + 5 * (1 - float((abs(data.loc[ikey]) - limits[0]) / (limits[1] - limits[0])))
                color = cmap(float((abs(data.loc[ikey]) - limits[0])/(limits[1] - limits[0])))

                ax.quiver(x_position, y_position, x_direction, y_direction, scale = scale-3, headwidth=4, headaxislength=2, headlength=2, pivot="mid", zorder=3, color="black", linewidth=0)
                # ax.quiver(x_position, y_position, x_direction, y_direction, scale = scale, headwidth=4, headaxislength=2, headlength=2, pivot="mid", zorder=4, color=color)

        else:
            print(f"Dataframe index {ikey} has no corresponding shape")


def plot_patches(ax, data, patch_coll, cmap, limits):
    for ikey in data.index.values:
        if ikey in patch_coll:
            ax.add_collection(
                PatchCollection(
                    patch_coll[ikey],
                    facecolor=cmap(float((data.loc[ikey] - limits[0]) / (limits[1] - limits[0]))),
                    edgecolor="k",
                    linewidths=0.3,
                    zorder=2,
                )
            )
        else:
            print(f"Dataframe index {ikey} has no corresponding shape")


def plot_background(ax, data, patch_coll):
    for ikey in data.index.values:
        if ikey in patch_coll:
            ax.add_collection(
                PatchCollection(
                    patch_coll[ikey],
                    facecolor="#fcfcfc",
                    edgecolor="k",
                    linewidths=0.3,
                    zorder=2,
                )
            )
        else:
            print(f"Dataframe index {ikey} has no corresponding shape")


def add_cbar(ax, limits, cmap, clabel):
    norm = matplotlib.colors.Normalize(vmin=limits[0], vmax=limits[1])
    scamap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

    cbar = plt.colorbar(scamap, ax=ax)

    if clabel is not None:
        cbar.set_label(clabel, rotation=90, labelpad=10)


def plot_choropleth(
    data, shp, shp_attr, ax=None, basemap=True, lat=None, lon=None, limits=None, cmap="viridis", cbar=True, clabel=None
):
    patch_coll = load_polypatches(shp, shp_attr)
    idx_match = [i for i in set(data.index.get_level_values(0)) if i in list(patch_coll.keys())]
    data = data.loc[idx[idx_match]]

    if ax is None:
        fig = plt.figure(figsize=[12, 8])
        ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])

    if basemap:
        add_basemap(lat, lon)

    cmap, limits = get_cmap(cmap, data, limits)

    plot_patches(ax, data, patch_coll, cmap, limits)

    if cbar:
        add_cbar(ax, limits, cmap, clabel)

    return fig, ax


def plot_network_caps(
    data,
    shp,
    shp_attr,
    ax=None,
    basemap=True,
    backgr_shp=None,
    backgr_attr=None,
    lat=None,
    lon=None,
    limits=None,
    cmap="viridis",
    cbar=True,
    clabel=None,
):
    patch_coll = load_polypatches(shp, shp_attr)
    if backgr_shp is not None and backgr_attr is not None:
        backgr_coll = load_polypatches(backgr_shp, backgr_attr)
    idx_match = [i for i in set(data.index.get_level_values(0)) if i in list(patch_coll.keys())]
    data = data.loc[idx[idx_match]]

    if ax is None:
        fig = plt.figure(figsize=[12, 8])
        ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])

    if basemap:
        add_basemap(lat, lon)

    cmap, limits = get_cmap(cmap, data, limits)

    if backgr_shp is not None and backgr_attr is not None:
        plot_background(ax, data, backgr_coll)

    plot_lines(ax, data, patch_coll, cmap, limits)

    if cbar:
        add_cbar(ax, limits, cmap, clabel)

    return fig, ax


def plot_network_flows(
    data,
    shp,
    shp_attr,
    ax=None,
    basemap=True,
    backgr_shp=None,
    backgr_attr=None,
    lat=None,
    lon=None,
    limits=None,
    cmap="viridis",
    cbar=True,
    clabel=None,
):
    patch_coll = load_polypatches(shp, shp_attr)
    if backgr_shp is not None and backgr_attr is not None:
        backgr_coll = load_polypatches(backgr_shp, backgr_attr)
    idx_match = [i for i in set(data.index.get_level_values(0)) if i in list(patch_coll.keys())]
    data = data.loc[idx[idx_match]]

    if ax is None:
        fig = plt.figure(figsize=[12, 8])
        ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])

    if basemap:
        add_basemap(lat, lon)

    cmap, limits = get_cmap(cmap, data, limits)

    if backgr_shp is not None and backgr_attr is not None:
        plot_background(ax, data, backgr_coll)

    plot_arrows(ax, data, patch_coll, cmap, limits)

    if cbar:
        add_cbar(ax, limits, cmap, clabel)

    return fig, ax


def filter_dataframe(df, filter):
    for idx, sel in filter.items():
        if idx in df.index.names:
            _indexer = [slice(None)] * len(df.index.levels)
            idx_pos = df.index.names.index(idx)
            _indexer[idx_pos] = sel
            indexer = tuple(_indexer)
            df = df.loc[indexer, slice(None)]
    return df


def generate_choropleths(cfg_choro, gdxeval=None):
    assert "gdx" in cfg_choro and "shp" in cfg_choro and "shp_attr" in cfg_choro
    shp = cfg_choro["shp"]
    shp_attr = cfg_choro["shp_attr"]
    lat = cfg_choro["lat"] if "lat" in cfg_choro else None
    lon = cfg_choro["lon"] if "lon" in cfg_choro else None

    if gdxeval is None:
        gdx = cfg_choro["gdx"]
        res_gdx = GDXEval(gdx)
    else:
        res_gdx = gdxeval

    if "plots" in cfg_choro:
        for k, cfg_plot in cfg_choro["plots"].items():
            assert "symbol" in cfg_plot
            scale = cfg_plot["scale"] if "scale" in cfg_plot else 1
            label = cfg_plot["label"] if "label" in cfg_plot else None

            data = res_gdx[cfg_plot["symbol"]] * float(scale)
            if "scenario" in data.index.names:
                data.index = data.index.droplevel("scenario")
            print(data)

            if "filter" in cfg_plot:
                data = filter_dataframe(data, cfg_plot["filter"])

            if "accNodesModel" in data.index.names:
                spatial = "accNodesModel"
            elif "nodesModel" in data.index.names:
                spatial = "nodesModel"

            if "accYears" in data.index.names:
                years = "accYears"
            elif "years" in data.index.names:
                years = "years"

            for y in set(data.index.get_level_values(years)):
                data_plt = filter_dataframe(data, {years: y})
                data_plt = data_plt.groupby(spatial).sum()
                print(f"Creating choropleth plot {k} for year {y}")
                fig, ax = plot_choropleth(data_plt, shp, shp_attr, lat=lat, lon=lon, clabel=label)
                plot_dir = Path(plotdir, "choropleth")
                plot_dir.mkdir(exist_ok=True, parents=True)
                plt.savefig(f"{plot_dir}/{k}_{y}.png", bbox_inches="tight", dpi=300)
                plt.savefig(f"{plot_dir}/{k}_{y}.svg", bbox_inches="tight")


def generate_networks(cfg_netw, gdxeval=None):
    assert "gdx" in cfg_netw and "shp" in cfg_netw and "shp_attr" in cfg_netw
    shp = cfg_netw["shp"]
    shp_attr = cfg_netw["shp_attr"]
    lat = cfg_netw["lat"] if "lat" in cfg_netw else None
    lon = cfg_netw["lon"] if "lon" in cfg_netw else None

    if gdxeval is None:
        gdx = cfg_netw["gdx"]
        res_gdx = GDXEval(gdx)
    else:
        res_gdx = gdxeval

    if "plots" in cfg_netw:
        for k, cfg_plot in cfg_netw["plots"].items():
            assert "symbol" in cfg_plot
            scale = cfg_plot["scale"] if "scale" in cfg_plot else 1
            label = cfg_plot["label"] if "label" in cfg_plot else None

            data = res_gdx[cfg_plot["symbol"]] * float(scale)
            if "scenario" in data.index.names:
                data.index = data.index.droplevel("scenario")
            print(data)

            data["link"] = (
                data.index.get_level_values(0).astype(str) + "__" + data.index.get_level_values(1).astype(str)
            )
            data = data.droplevel(["nodesModel_start", "nodesModel_end", "linksModel"])
            data = data.set_index("link", append=True)

            if "filter" in cfg_plot:
                data = filter_dataframe(data, cfg_plot["filter"])

            if "accYears" in data.index.names:
                years = "accYears"
            elif "years" in data.index.names:
                years = "years"

            for y in set(data.index.get_level_values(years)):
                data_plt = filter_dataframe(data, {years: y})
                data_plt = data_plt.groupby("link").sum().abs()

                print(f"Creating network plot {k} for year {y}")
                fig, ax = plot_network_caps(
                    data_plt, shp, shp_attr, lat=lat, lon=lon
                )  # backgr_shp=backgr_shp, backgr_attr=backgr_attr,
                plot_dir = Path(plotdir, "network")
                plot_dir.mkdir(exist_ok=True, parents=True)
                plt.savefig(f"{plot_dir}/{k}_{y}.png", bbox_inches="tight", dpi=300)
                plt.savefig(f"{plot_dir}/{k}_{y}.svg", bbox_inches="tight")


def run_plotting(yml, gdxeval=None):
    with open(yml, "r") as yf:
        cfg = yaml.load(yf, Loader=yaml.FullLoader)

    if "choropleth" in cfg:
        generate_choropleths(cfg["choropleth"], gdxeval)

    if "network" in cfg:
        generate_networks(cfg["network"], gdxeval)


if __name__ == "__main__":
    # gdx = "../results/myopic-debug/remix.gdx"
    # plotting defaults and options
    plotdir = os.path.abspath("plots_H2") + "/"

    # gdx = "../results/myopic-debug/co2_targets/20Gt/remix.gdx"
    gdx = {"2030": "../results/H2/co2_targets/20Gt/targets/2030/remix.gdx",
           "2040": "../results/H2/co2_targets/20Gt/targets/2040/remix.gdx",
           "2045": "../results/H2/co2_targets/20Gt/targets/2045/remix.gdx",
           "2050": "../results/H2/co2_targets/20Gt/targets/2050/remix.gdx"}

    rename_years = {i: str(i) for i in range(1900, 2100)}

    rename_techs = {
        "CCGT": "CCGT",
        "CCGT_H2": "CCGT (H2)",
        "csp_powerblock": "CSP",
        "pv_central_fixed": "PV open area (fixed)",
        "pv_central_track_azimuth": "PV open area (tracking)",
        "wind_onshore": "Onshore wind",
        "wind_offshore_floating": "Offshore wind (floating)",
        "wind_offshore_foundation": "Offshore wind (foundation)",
    }

    results = GDXEval(gdx)
    run_plotting("nagsys_path.yaml", results)

    cba = results["commodity_balance_annual"].rename(index=rename_years)
    cba = cba[cba.abs() > 0.01].dropna()
    if "scenario" in cba.index.names:
        cba.index = cba.index.droplevel("scenario")

    cap = results["converter_caps"].rename(index=rename_years)
    cap = cap[cap.abs() > 0.01].dropna()
    if "scenario" in cap.index.names:
        cap.index = cap.index.droplevel("scenario")

    ind_dtl = results["indicator_accounting_detailed"].rename(index=rename_years)
    if "scenario" in ind_dtl.index.names:
        ind_dtl.index = ind_dtl.index.droplevel("scenario")

    techs = list(set(cba.index.get_level_values(2)))
    techs_pv = [i for i in techs if i.lower().startswith("pv")]
    techs_wind = [i for i in techs if i.lower().startswith("wind")]
    nodes = [n for n in cba.index.get_level_values(0) if n != "global"]
    rename_nodes = {n: "_".join(n.split("_")[:-1]) for n in nodes}

    elec = cba.loc[idx[nodes, :, :, "Elec", "netto"]]
    ch4 = cba.loc[idx[nodes, :, :, "CH4", "netto"]]
    h2 = cba.loc[idx[nodes, :, :, "H2", "netto"]]


    # print(ind.loc[idx["CO2Emission", :, :, :]].div(1e3).sort_values("value"))
    # print(cba.loc[idx[:,:,:,"Elec","netto"]].div(1e3).unstack("accNodesModel").sort_values(("value", "global")))
    # print(cba.loc[idx[:,:,:,"H2","netto"]].div(1e3).unstack("accNodesModel").sort_values(("value", "global")))
    # print(cba.loc[idx[:,:,:,"CH4","netto"]].div(1e3).unstack("accNodesModel").sort_values(("value", "global")))
    # print(cba.loc[idx[:,:,:,"LNG","netto"]].div(1e3).unstack("accNodesModel").sort_values(("value", "global")))
    # print(cba.loc[idx[:,:,:,"LH2","netto"]].div(1e3).unstack("accNodesModel").sort_values(("value", "global")))
    # print(cba.loc[idx[:,:,:,"Elec_battery","brutto"]].unstack("accNodesModel").sort_values(("value", "global")).T.dropna())
    # print(cba.loc[idx["global",:,["Households", "Industry", "Tertiary", "Transport"],:,"netto"]].div(1e3).sort_values("value"))

    # print(cap.loc[idx[:,:,:,"Elec","total"]].unstack("accNodesModel").sort_values(("value", "global")))
    # print(cap.loc[idx[:,:,:,"H2","total"]].unstack("accNodesModel").sort_values(("value", "global")))
    # print(cap.loc[idx[:,:,:,"CH4","total"]].unstack("accNodesModel").sort_values(("value", "global")))

    # print(cba.loc[idx[:, :, ["csp_powerblock"], ["Elec"], ["netto"]]].sort_values("value"))

    # # Plotting test
    # elec = elec.groupby(["accYears", "techs"]).sum().div(1e3)
    # # elec = elec[elec > 0.1].dropna().fillna(0)
    # # elec = elec[elec < 0.1].dropna().fillna(0)
    # fig = px.bar(elec.reset_index(), x="accYears", y="value", color="techs")
    # fig.show()

    # elec = cba.loc[idx[nodes, :, :, "H2", "netto"]]
    # elec = elec.groupby(["accYears", "techs"]).sum().div(1e3)
    # # elec = elec[elec > 0.1].dropna().fillna(0)
    # # elec = elec[elec < 0.1].dropna().fillna(0)
    # fig = px.bar(elec.reset_index(), x="accYears", y="value", color="techs")
    # fig.show()

    # elec = cba.loc[idx[nodes, :, :, "CH4", "netto"]]
    # elec = elec.groupby(["accYears", "techs"]).sum().div(1e3)
    # # elec = elec[elec > 0.1].dropna().fillna(0)
    # # elec = elec[elec < 0.1].dropna().fillna(0)
    # fig = px.bar(elec.reset_index(), x="accYears", y="value", color="techs")
    # fig.show()

    # elec = cap.loc[idx[nodes, :, :, "Elec", "total"]]
    # elec = elec.groupby(["accYears", "techs"]).sum()
    # fig = px.bar(elec.reset_index(), x="accYears", y="value", color="techs")
    # fig.show()

    # co2 = ind_dtl.loc[idx["CO2Emission", :, :, :]]
    # co2 = co2.groupby(["years", "techs"]).sum()
    # fig = px.bar(co2.reset_index(), x="years", y="value", color="techs")
    # fig.show()


    # gen_per_tech = gen.loc[idx[nodes_model,:,:]].rename(index=rename_techs).groupby(["accYears", "techs"]).sum().unstack(["techs"])
    # gen_per_tech.columns = gen_per_tech.columns.get_level_values(1)
    # fig = plt.figure(figsize=[8, 6])
    # ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])
    # gen_per_tech.plot(ax=ax, kind='bar', stacked=True)
    # plt.legend(loc=(1.05,0.1))
    # ax.set_xlabel("Year")
    # ax.set_ylabel("Annual power generation in TWh")
    # plt.xticks(rotation=0)
    # plt.savefig(f"{plotdir}/electricity_gen_per_tech.png", bbox_inches='tight', dpi=300)

    elec = elec.groupby(["accYears", "techs"]).sum().div(1e3)

    demand_per_tech = elec.rename(index=rename_techs).unstack(["techs"])
    demand_per_tech.columns = demand_per_tech.columns.get_level_values(1)
    fig = plt.figure(figsize=[5, 7])
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])
    demand_per_tech.plot(ax=ax, kind='bar', stacked=True)
    plt.legend(loc=(1.05,0.0))
    ax.set_xlabel("Year")
    ax.set_ylabel("Annual electricity balance in TWh")
    plt.xticks(rotation=0)
    plt.savefig(f"{plotdir}/electricity_balance.png", bbox_inches='tight', dpi=300)


    ch4 = ch4.groupby(["accYears", "techs"]).sum().div(1e3)
    demand_per_tech = ch4.rename(index=rename_techs).unstack(["techs"])
    demand_per_tech.columns = demand_per_tech.columns.get_level_values(1)
    fig = plt.figure(figsize=[5, 7])
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])
    demand_per_tech.plot(ax=ax, kind='bar', stacked=True)
    plt.legend(loc=(1.05,0.05))
    ax.set_xlabel("Year")
    ax.set_ylabel("Annual methane balance in TWh")
    plt.xticks(rotation=0)
    plt.savefig(f"{plotdir}/methane_balance.png", bbox_inches='tight', dpi=300)



    h2 = h2.groupby(["accYears", "techs"]).sum().div(1e3)
    demand_per_tech = h2.rename(index=rename_techs).unstack(["techs"])
    demand_per_tech.columns = demand_per_tech.columns.get_level_values(1)
    fig = plt.figure(figsize=[5, 7])
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])
    demand_per_tech.plot(ax=ax, kind='bar', stacked=True)
    plt.legend(loc=(1.05,0.05))
    ax.set_xlabel("Year")
    ax.set_ylabel("Annual hydrogen balance in TWh")
    plt.xticks(rotation=0)
    plt.savefig(f"{plotdir}/hydrogen_balance.png", bbox_inches='tight', dpi=300)

    # print(h2.groupby(["accYears", "techs"]).sum().div(1e3))