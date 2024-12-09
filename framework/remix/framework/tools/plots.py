import yaml
import fiona
import pyproj
import numpy as np
import pandas as pd
from pathlib import Path
from shapely.validation import make_valid
from shapely.geometry import Polygon, MultiPolygon, LineString, shape
from shapely.ops import transform
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon as PolyPatch
from matplotlib.collections import PatchCollection, LineCollection
from remix.framework.tools.gdx import GDXEval
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import itertools

from typing import Optional

idx = pd.IndexSlice


def make_geo_valid(ob):
    return ob if ob.is_valid else make_valid(ob.geoms)


def load_polypatches(shp, shp_attr, crs=3857):
    coll = {}
    with fiona.open(shp, "r") as f:
        project = transform_from_latlon(crs)
        for i in f:
            geo = shape(i["geometry"])
            geo_t = transform(project, geo)
            coll[i["properties"][shp_attr]] = geo_t

    patch_coll = {}
    for k, i in coll.items():
        if isinstance(i, MultiPolygon):
            patch_coll[k] = [PolyPatch(xy=np.array(j.exterior.coords.xy).T, closed=True) for j in i.geoms]
        elif isinstance(i, Polygon):
            patch_coll[k] = [PolyPatch(xy=np.array(i.exterior.coords.xy).T, closed=True)]
        elif isinstance(i, LineString):
            patch_coll[k] = [np.array(i.coords.xy).T]
            if "__" in k:
                k_r = "__".join(reversed(k.split("__")))
                patch_coll[k_r] = [np.fliplr(np.array(i.coords.xy)).T]

    return patch_coll


def get_cmap(cmap, data, limits: Optional[list[int | float]] = None):
    cmap = matplotlib.cm.get_cmap(cmap)
    if limits is None:
        limits = [None, None]
    if limits[0] is None:
        limits[0] = data.min().iloc[0]
    if limits[1] is None:
        limits[1] = data.max().iloc[0]
    return cmap, limits


def build_linestrings(centroids: dict):
    return {(i, j): LineString((i_pos, j_pos)) for i, i_pos in centroids.items() for j, j_pos in centroids.items() if i != j}


def transform_from_latlon(crs: int = 3857):
    crs_source = pyproj.CRS("EPSG:4326")  # lat-lon coord system
    crs_target = pyproj.CRS(f"EPSG:{crs}")

    return pyproj.Transformer.from_crs(crs_source, crs_target, always_xy=True).transform


def plot_lines(ax, data: pd.DataFrame, patch_coll, cmap, limits: Optional[list[int | float]] = None, crs: int = 3857):
    trans = transform_from_latlon(crs)
    for ikey in data.index.values:
        if ikey in patch_coll:
            if limits is not None:
                edgecolor = cmap((data.loc[ikey].iloc[0] - limits[0]) / (limits[1] - limits[0]))
                lwidth = (data.loc[ikey].iloc[0] - limits[0]) / (limits[1] - limits[0])
            else:
                edgecolor = "Reds"
                lwidth = 1/8
            seg_in_metres = transform(trans, patch_coll[ikey])
            line_segments = LineCollection(
                    segments=[list(seg_in_metres.coords)],
                    edgecolor=edgecolor,
                    linewidths=8 * lwidth,
                    zorder=2,
                )

            ax.add_collection(line_segments)
        else:
            print(f"DataFrame index {ikey} has no corresponding shape")

    return ax


def plot_patches(ax, data, patch_coll, cmap, limits: Optional[list[int | float]] = None):
    for ikey in data.index.values:
        if ikey in patch_coll:
            if limits is not None:
                facecolor = cmap((data.loc[ikey, "value"] - limits[0]) / (limits[1] - limits[0]))
            else:
                facecolor = "Reds"

            ax.add_collection(
                PatchCollection(
                    patch_coll[ikey],
                    facecolor=facecolor,
                    edgecolor="k",
                    linewidths=0.3,
                    zorder=2,
                )
            )
        else:
            print(f"DataFrame index {ikey} has no corresponding shape")

    return ax


def plot_background(ax, data, patch_coll):
    nodes = set(itertools.chain.from_iterable(data.index))
    for ikey in nodes:
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
            print(f"DataFrame index {ikey} has no corresponding shape")


def plot_choropleth(
    data,
    shp,
    shp_attr,
    fig=None,
    ax=None,
    title: Optional[str] = None,
    background: bool = True,
    lat=None,
    lon=None,
    limits: Optional[list[int | float]] = None,
    cmap: str = "viridis",
    cbar: bool = True,
    clabel: Optional[str] = None,
    patch_coll=None,
    crs: int = 3857,
):
    if patch_coll is None:
        patch_coll = load_polypatches(shp, shp_attr, crs)

    idx_match = [i for i in set(data.index.get_level_values(0)) if i in list(patch_coll.keys())]
    data = data.loc[idx[idx_match]]

    if ax is None:
        fig = plt.figure(figsize=[6, 8])
        ax = fig.add_axes([0.02, 0.02, 0.96, 0.96], projection=ccrs.epsg(crs))

    if background:
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor="#e0e0e0"))
        # ax.add_feature(cfeature.NaturalEarthFeature('physical', 'lakes', '10m', edgecolor='face', facecolor="#ffffff"))

    cmap, limits = get_cmap(cmap, data, limits)

    if title is not None:
        ax.set_title(title)

    ax = plot_patches(ax, data, patch_coll, cmap, limits)

    if cbar:
        cax1 = ax.inset_axes([1.03, 0.0, 0.07, 1.0])
        norm = matplotlib.colors.Normalize(vmin=limits[0], vmax=limits[1]) if limits is not None else None
        scamap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        cbar = plt.colorbar(scamap, cax=cax1, orientation="vertical", label=clabel)

    ax.set_extent([*lon, *lat])

    return fig, ax


def plot_network(
    data,
    shp=None,
    shp_attr=None,
    centroids: Optional[dict] = None,
    fig=None,
    ax=None,
    title: Optional[str] = None,
    background: bool = True,
    backgr_shp=None,
    backgr_attr=None,
    lat=None,
    lon=None,
    limits: Optional[list[int | float]] = None,
    cmap: str = "viridis",
    cbar=True,
    clabel=None,
    patch_coll=None,
    crs: int = 3857,
):
    # Take abs for undirected flow volumes
    data = data.abs()

    if patch_coll is None:
        if centroids is not None:
            patch_coll = build_linestrings(centroids)

        if (shp is not None) and (shp_attr is not None):
            patch_coll = load_polypatches(shp, shp_attr, crs)

            if not any("__" in s for s in patch_coll):
                patch_coll = build_linestrings(centroids)

        if centroids is None and shp is None and shp_attr is None:
            raise "Need either shapefile or a location dictionary to generate network map."

    if backgr_shp is not None and backgr_attr is not None:
        backgr_coll = load_polypatches(backgr_shp, backgr_attr, crs)

    idx_match = [i for i in set(data.index.to_flat_index()) if i in list(patch_coll.keys())]
    data = data.loc[idx[idx_match]]

    if ax is None:
        fig = plt.figure(figsize=[6, 8])
        ax = fig.add_axes([0.02, 0.02, 0.96, 0.96], projection=ccrs.epsg(crs))

    if background:
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor="#e0e0e0"))
        # ax.add_feature(cfeature.NaturalEarthFeature('physical', 'lakes', '10m', edgecolor='face', facecolor="#ffffff"))

    cmap, limits = get_cmap(cmap, data, limits)

    if backgr_shp is not None and backgr_attr is not None:
        plot_background(ax, data, backgr_coll)

    if title is not None:
        ax.set_title(title)

    ax = plot_lines(ax, data.fillna(0), patch_coll, cmap, limits, crs)

    if cbar:
        cax1 = ax.inset_axes([1.03, 0.0, 0.07, 1.0])
        norm = matplotlib.colors.Normalize(vmin=limits[0], vmax=limits[1]) if limits is not None else None
        scamap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        cbar = plt.colorbar(scamap, cax=cax1, orientation="vertical", label=clabel)

    ax.set_extent([*lon, *lat])

    return fig, ax


def filter_dataframe(df: pd.DataFrame, filter) -> pd.DataFrame:
    for idx, sel in filter.items():
        if idx in df.index.names:
            _indexer = [slice(None)] * len(df.index.levels)
            idx_pos = df.index.names.index(idx)
            _indexer[idx_pos] = sel
            indexer = tuple(_indexer)
            df = df.loc[indexer, slice(None)]
    return df


def generate_choropleths(cfg_choro: dict, gdxeval=None, crs: int = 3857) -> None:
    assert "gdx" in cfg_choro and "shp" in cfg_choro and "shp_attr" in cfg_choro
    gdx = cfg_choro["gdx"]
    shp = cfg_choro["shp"]
    shp_attr = cfg_choro["shp_attr"]
    patch_coll = load_polypatches(shp, shp_attr, crs)

    if gdxeval is None:
        res_gdx = GDXEval(gdx)
    else:
        res_gdx = gdxeval

    lat = cfg_choro["lat"] if "lat" in cfg_choro else None
    lon = cfg_choro["lon"] if "lon" in cfg_choro else None

    if "plots" in cfg_choro and cfg_choro["plots"] is not None:
        for k, cfg_plot in cfg_choro["plots"].items():
            assert "symbol" in cfg_plot
            scale = cfg_plot["scale"] if "scale" in cfg_plot else 1
            label = cfg_plot["label"] if "label" in cfg_plot else None

            data = res_gdx[cfg_plot["symbol"]] * float(scale)

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

                # print(f"Creating choropleth plot {k} for year {y}")
                fig, ax = plot_choropleth(data_plt, shp, shp_attr, lat=lat, lon=lon, clabel=label, patch_coll=patch_coll)
                plot_dir = Path("plots", "choropleth")
                plot_dir.mkdir(exist_ok=True, parents=True)
                plt.savefig(f"{plot_dir}/{k}_{y}.png", bbox_inches="tight", dpi=300)
                plt.savefig(f"{plot_dir}/{k}_{y}.svg", bbox_inches="tight")
                print(f"Saved choropleth plot {plot_dir}/{k}_{y}.png")


def generate_networks(cfg_netw: dict, gdxeval=None, crs: int = 3857) -> None:
    assert "gdx" in cfg_netw and "shp" in cfg_netw and "shp_attr" in cfg_netw
    gdx = cfg_netw["gdx"]
    shp = cfg_netw["shp"]
    shp_attr = cfg_netw["shp_attr"]
    patch_coll = load_polypatches(shp, shp_attr, crs)

    if gdxeval is None:
        res_gdx = GDXEval(gdx)
    else:
        res_gdx = gdxeval

    lat = cfg_netw["lat"] if "lat" in cfg_netw else None
    lon = cfg_netw["lon"] if "lon" in cfg_netw else None

    if "plots" in cfg_netw:
        for k, cfg_plot in cfg_netw["plots"].items():
            assert "symbol" in cfg_plot
            scale = cfg_plot["scale"] if "scale" in cfg_plot else 1
            label = cfg_plot["label"] if "label" in cfg_plot else None

            data = res_gdx[cfg_plot["symbol"]] * float(scale)
            data["link"] = data.index.get_level_values(0).astype(str) + "__" + data.index.get_level_values(1).astype(str)
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
                data_plt = data_plt.groupby("link").sum()

                # print(f"Creating network plot {k} for year {y}")
                fig, ax = plot_network(data_plt, shp, shp_attr, lat=lat, lon=lon, patch_coll=patch_coll)#, backgr_shp=backgr_shp, backgr_attr=backgr_attr)
                plot_dir = Path("plots", "network")
                plot_dir.mkdir(exist_ok=True, parents=True)
                plt.savefig("test.png", bbox_inches="tight", dpi=300)
                plt.savefig(f"{plot_dir}/{k}_{y}.png", bbox_inches="tight", dpi=300)
                plt.savefig(f"{plot_dir}/{k}_{y}.svg", bbox_inches="tight")
                print(f"Saved network plot {plot_dir}/{k}_{y}.png")


def run_plotting(yml: str) -> None:
    with open(yml, "r") as yf:
        cfg = yaml.load(yf, Loader=yaml.FullLoader)

    if "choropleth" in cfg:
        generate_choropleths(cfg["choropleth"])

    if "network" in cfg:
        generate_networks(cfg["network"])
