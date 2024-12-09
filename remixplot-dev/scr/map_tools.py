import math
import numpy as np
import pandas as pd
import geopandas as gpd
import plotly.graph_objects as go

from dataclasses import dataclass
from shapely.geometry import Polygon, MultiPolygon


@dataclass
class Fig:

    @staticmethod
    def fig_layout(zoom: int) -> go.Figure:
        fig = go.Figure()
        fig.update(
            layout=go.Layout(
                showlegend=False,
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                margin={"l": 0, "t": 35, "b": 0, "r": 0},
            ),
        ).update_layout(
            geo=dict(
                projection_type="eckert4",
                projection_scale=zoom,
                resolution=50,
                showland=True,
                landcolor="white",
                showocean=True,
                oceancolor="azure",
                showlakes=True,
                lakecolor="azure",
                showcountries=True,
                countrycolor="light gray",
                showcoastlines=True,
                coastlinecolor="black",
            )
        )

        return fig

    @staticmethod
    def fig_country(
        lat: list[float],
        lon: list[float],
        node_name: str,
        color: str | dict = "gray",
        key: str | None = None,
    ) -> go.Scattergeo:

        if isinstance(color, dict):
            color = color[key]
        return go.Scattergeo(
            lon=lon,
            lat=lat,
            mode="lines",
            line=dict(width=2, color=color),
            fill="toself",
            fillcolor=color,
            opacity=0.5,
            name=node_name,
        )

    @staticmethod
    def fig_nodes(node_lat: float, node_lon: float, node_name: str) -> go.Scattergeo:

        return go.Scattergeo(
            mode="markers",
            lon=[node_lon],
            lat=[node_lat],
            name=node_name,
            marker={"size": 10, "color": "black"},
        )

    @staticmethod
    def fig_trans_connections(
        node_dict: dict,
        node_start: str,
        node_end: str,
        ratioVal: float,
        width_fac: float,
    ) -> go.Scattergeo:

        lon_start, lat_start, lon_end, lat_end = lats_lons_nodes(
            node_dict, node_start, node_end
        )

        conn_width: float = max(1, ratioVal * width_fac)
        conn_color: str = "blue"

        return go.Scattergeo(
            lon=[lon_start, lon_end],
            lat=[lat_start, lat_end],
            mode="lines",
            line={"color": conn_color, "width": conn_width},
        )

    @staticmethod
    def fig_trans_arrows(
        node_dict: dict,
        node_start: str,
        node_end: str,
        ratioVal: float,
        val: float,
        width_fac: float,
    ) -> go.Scattergeo:

        arrow_fac: float = 3

        lon_start, lat_start, lon_end, lat_end = lats_lons_nodes(
            node_dict, node_start, node_end
        )

        lat_symbol, lon_symbol = get_symbol_distances(
            lat_start, lat_end, lon_start, lon_end
        )
        symbol_size: float = max(2 * arrow_fac, ratioVal * width_fac * arrow_fac)
        symbol_theta = get_angle(lat_start, lat_end, lon_start, lon_end, val)

        conn_color: str = "azure"
        return go.Scattergeo(
            lon=lon_symbol,
            lat=lat_symbol,
            mode="markers",
            marker=dict(
                symbol="arrow", size=symbol_size, angle=symbol_theta, color=conn_color
            ),
        )

    @staticmethod
    def fig_connections(
        node_dict: dict, node_start: str, node_end: str
    ) -> go.Scattergeo:

        lon_start, lat_start, lon_end, lat_end = lats_lons_nodes(
            node_dict, node_start, node_end
        )

        conn_width = 1.5
        conn_color = "Red"
        return go.Scattergeo(
            lon=[lon_start, lon_end],
            lat=[lat_start, lat_end],
            mode="lines",
            line={"color": conn_color, "width": conn_width},
        )


def get_map_center(fig: go.Figure, node_dict: dict) -> go.Figure:

    return fig.update_layout(
        geo=dict(
            center=dict(
                lon=np.mean([val[0] for val in node_dict.values()]),
                lat=np.mean([val[1] for val in node_dict.values()]),
            )
        )
    )


def get_contour_coordinates(df_filt: gpd.GeoDataFrame) -> tuple[list[str], list[str]]:

    all_lat_coords, all_lon_coords = [], []

    for _, row in df_filt.iterrows():
        geometry = row["geometry"]

        if isinstance(geometry, Polygon):
            lat, lon = geometry.exterior.xy
            all_lat_coords.extend(lat)
            all_lon_coords.extend(lon)

        elif isinstance(geometry, MultiPolygon):
            for poly in geometry.geoms:
                lat, lon = poly.exterior.xy
                all_lat_coords.extend(lat)
                all_lon_coords.extend(lon)
                all_lat_coords.append(None)
                all_lon_coords.append(None)

    return all_lon_coords, all_lat_coords


def calculate_centroid(lats: list[float], lons: list[float]) -> tuple[float, float]:

    lats = list(filter(None, lats))
    lons = list(filter(None, lons))
    lat = sum(lats) / len(lats)
    lon = sum(lons) / len(lons)

    return (lat, lon)


def get_center(
    lats: list[float], lons: list[float], id: str
) -> tuple[list[float], list[float]]:

    if id == "ES":
        lon_cent, lat_cent = -3.703533, 40.417065
    elif id == "FR":
        lon_cent, lat_cent = 2.2639, 46.4609
    elif id == "PT":
        lon_cent, lat_cent = -8.145690, 39.675621
    elif id == "DE6":
        lon_cent, lat_cent = 9.989791, 53.553222
    elif id == "SE":
        lon_cent, lat_cent = 14.612639, 59.707775
    elif id == "NO":
        lon_cent, lat_cent = 8.639684, 60.899451
    elif id == "IT":
        lon_cent, lat_cent = 11.069551, 44.451735
    elif id == "DK_E":
        lon_cent, lat_cent = 11.762922, 55.363165
    else:
        lat_cent, lon_cent = calculate_centroid(lats, lons)

    return lat_cent, lon_cent


def get_symbol_distances(
    lat_start: float, lat_end: float, lon_start: float, lon_end: float
) -> tuple[list[float], list[float]]:

    lon_scatter, lat_scatter = [], []
    loni, lati = lon_start, lat_start

    diff_lat = lat_end - lat_start
    diff_lon = lon_end - lon_start

    flowDir_lat = 1 if lat_start < lat_end else -1
    flowDir_lon = 1 if lon_start < lon_end else -1

    for _ in range(3):
        lati += flowDir_lat * abs(diff_lat) / 4
        loni += flowDir_lon * abs(diff_lon) / 4

        lat_scatter += [lati]
        lon_scatter += [loni]

    return lat_scatter, lon_scatter


def get_angle(
    lat_start: float, lat_end: float, lon_start: float, lon_end: float, val: float
) -> float:

    theta0: float = math.atan2(lon_end - lon_start, lat_end - lat_start)
    theta: float = math.degrees(theta0)
    theta: float = theta if val > 0 else theta + 180

    return theta


def lats_lons_nodes(
    node_dict: dict, node_start: str, node_end: str
) -> tuple[float, float, float, float]:

    lon_start, lat_start = node_dict[node_start]
    lon_end, lat_end = node_dict[node_end]

    return lon_start, lat_start, lon_end, lat_end
