import random
import matplotlib.colors as mpl_colors


class Color_schema:
    def __init__(self, utils) -> None:

        self.blue: str = "#3944BC"
        self.green: str = "#98BF64"
        self.dark_gray: str = "#232023"

        self.defined_colors: dict = {
            # -- wind --
            "Wind": "#9dc3e6",
            "Wind_off": "#9dc3e6",
            "WindOffshore": "#9dc3e6",
            "Wind_on": "#9dc3e6",
            "WindOnshore": "#9dc3e6",
            # -- solar --
            "Solar": "#ffc000",
            "PV": "#ffc000",
            "Pv": "#ffc000",
            "Hydro": "#0070c0",
            # -- elec --
            "Electricity": "#70ad47",
            "Elec": "#70ad47",
            "Stored_LiIon": "#cc99ff",
            "Elec_stored": "#cc99ff",
            # -- heat --
            "Heat": "#b40000",
            # -- co2 --
            "CO2": "#76583a",
            "Co2": "#76583a",
            # -- fuels --
            "CH4": "#ff5f2d",
            "Coal": "#8c4471",
            "Lignite": "#b67878",
            "Nuclear": "#651278",
            "Oil": "#4c0000",
            # -- renew fuel --
            "Geothermal": "#a9d18e",
            "Biomass": "#ffd279",
            "Biogas": "#e7864b",
            "Hydrogen": "#439b80",
            "Pumped_Hydro": "#002060",
            # --
            "Waste": "#864300",
            "Compressed_air": "#33cccc",
        }

        self.tech_color_dict: dict = self.get_tech_color_dict(utils)

    def get_tech_color_dict(self, utils) -> dict:

        colors_dict: dict = dict(mpl_colors.cnames.items())
        hex_colors: list[str] = list(colors_dict.values())

        tech_color_dict: dict = {
            tec: random.choice(hex_colors)
            for tec in (utils.all_techs + utils.all_tech_groups)
        }

        for tec in tech_color_dict.keys():
            if tec in self.defined_colors.keys():
                tech_color_dict[tec] = self.defined_colors[tec]

        return tech_color_dict
