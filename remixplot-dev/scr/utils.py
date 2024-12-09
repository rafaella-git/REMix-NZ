import json
import pandas as pd
import networkx as nx
import geopandas as gpd

from remix.framework.tools.gdx import GDXEval


class Dataprep_utils:

    def __init__(
        self,
        gdx_path: str,
        geojson_path: str | None = None,
        mapping_path: str | None = None,
    ) -> None:

        self.add_dataframes(gdx_path)
        self.add_graph(gdx_path)
        self.add_geojson_mapping(geojson_path, mapping_path)

    def add_dataframes(self, gdx_path: str) -> None:
        gdx: GDXEval = GDXEval(gdx_path, fill_ts=False)

        # -- indicator ------------------------------------------------------------------
        self.indicator_accounting_comp: pd.DataFrame = gdx["indicator_accounting_comp"]

        # -- techs ----------------------------------------------------------------------
        self.all_techs: list[str] = list(gdx["techs"].reset_index()["uni"])
        self.all_tech_groups: list[str] = list(
            set([tec.split("_")[0] for tec in self.all_techs])
        )

        # -- timemodel ------------------------------------------------------------------
        timemodel_df: pd.DataFrame = gdx["timeModeltoCalc"]
        timemodel_df.index = timemodel_df.index.str.replace("tm", "").astype(int)
        self.timemodel_df: pd.DataFrame = timemodel_df

        # -- commodity_balance ----------------------------------------------------------
        self.commodity_balance_df: pd.DataFrame = gdx["commodity_balance"]

        if not self.commodity_balance_df.empty:
            self.generation_df: pd.DataFrame = self.commodity_balance_df[
                self.commodity_balance_df["value"] > 0
            ]
            self.demand_df: pd.DataFrame = self.commodity_balance_df[
                self.commodity_balance_df["value"] < 0
            ]

        # -- annual commodity_balance ---------------------------------------------------
        self.commodity_balance_annual_df: pd.DataFrame = gdx["commodity_balance_annual"]

        if not self.commodity_balance_annual_df.empty:
            self.generation_annual_df: pd.DataFrame = self.commodity_balance_annual_df[
                self.commodity_balance_annual_df["value"] > 0
            ]
            self.demand_annual_df: pd.DataFrame = self.commodity_balance_annual_df[
                self.commodity_balance_annual_df["value"] < 0
            ]

        # -- converter ------------------------------------------------------------------
        self.converter_units_df: pd.DataFrame = gdx["converter_units"]
        self.converter_caps_df: pd.DataFrame = gdx["converter_caps"]

        # -- storage --------------------------------------------------------------------
        self.storage_units_df: pd.DataFrame = (
            gdx["storage_units"].droplevel(3)
            if "vintage" in gdx["storage_units"].index.names
            else gdx["storage_units"]
        )
        self.storage_level_df: pd.DataFrame = gdx["storage_level_out"]

        # -- transfer_flow --------------------------------------------------------------
        self.trans_flow_df: pd.DataFrame = gdx["transfer_flows_annual"].reset_index()

    def add_graph(self, gdx_path: str) -> None:
        gdx = GDXEval(gdx_path)

        try:
            links_df = gdx["map_linksModel"].reset_index()
            links_df = links_df.drop("element_text", axis=1)
            links_data: list = links_df["linksData"]
            self.graph = nx.Graph([st_end.split("__") for st_end in links_data])
        except:
            pass

    def add_geojson_mapping(
        self, geojson_path: str | None, mapping_path: str | None
    ) -> None:

        if geojson_path is not None:
            self.geojson: gpd.GeoDataFrame = gpd.read_file(geojson_path)

        if mapping_path is not None:
            self.mapping: dict = json.load(open(mapping_path))
        else:
            try:
                self.mapping = {node: node for node in self.graph.nodes()}
            except:
                pass
