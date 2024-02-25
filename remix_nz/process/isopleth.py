import pandas as pd
import seaborn as sns
import matplotlib as plt
from remix.framework.tools.gdx import GDXEval
from remix.framework.tools.plots import plot_network, plot_choropleth
from pathlib import Path

sns.set(color_codes=True)
sns.set_style("whitegrid")
idx = pd.IndexSlice

name="nz-elec_2050"
path_base = "C:/Local/REMix"  
path_input = f"{path_base}/remix_nz/input"
path_output = f"{path_base}/remix_nz/output" 
path_result = f"{path_output}/{name}/result" 
path_profiles = f"{path_input}/profiles"
res = GDXEval(f"{path_result}/{name}.gdx")

print("hi")


def carnot_isopleths_gen(res):
    # Check if the path_result directory exists
    print(f"Path to result directory: {path_result}")
    if not Path(path_result).exists():
        print("Error: Result directory does not exist")

    # Check if the re_inst_csv file exists
    print(f"Path to profile data: {Path(path_profiles).joinpath('region_statistics_2012.csv')}")
    if not Path(path_profiles).joinpath('region_statistics_2012.csv').exists():
        print("Error: Profile data file does not exist")

    re_inst_csv = pd.read_csv(Path(path_profiles).joinpath("region_statistics_2012.csv"), index_col=[0, 1])
    re_techs = list(set(re_inst_csv.index.get_level_values(1)))
    pv_techs = [i for i in re_techs if i.startswith("pv")]
    csp_techs = []#[i for i in re_techs if i.startswith("csp")]
    wind_techs = [i for i in re_techs if i.startswith("wind")]

    map_dict = {
                # "Electrolyzer": {"title":"Electrolyzer output in GW", "techs":"Electrolyzer", "commodity":"H2", "scaling":1, "cmap": "viridis"},
    #             "Methaniser": {"title":"Methanizer output in GW", "techs":"Methaniser", "commodity":"CH4", "scaling":1, "cmap": "viridis"},
                "PV": {"title":"PV generation in GW", "techs":pv_techs, "commodity":"Elec", "scaling":1, "cmap": "viridis", "text": "b)"},
                "Wind": {"title":"Wind generation in GW", "techs":wind_techs, "commodity":"Elec", "scaling":1, "cmap": "viridis", "text": "a)"},
                # "CSP": {"title":"CSP generation in GW", "techs":"CSP_Powerblock", "commodity":"Elec", "scaling":1, "cmap": "viridis"},
    #            "RES": {"title":"Renewable generation in GW", "techs":["Photovoltaic","WindOnshore","WindOffshore","CSP_Powerblock"], "commodity":"Elec", "scaling":1, "cmap": "viridis"},
    #            "CHP": {"title":"CHP generation in GW", "techs":['TH_ExCCGT_NGas_XL','DH_Engine_NGas_M','DH_Engine_NGas_M', 'DH_FuelCell_H2_M'], "commodity":"Elec", "scaling":1, "cmap": "viridis"},
    #            "GT": {"title":"Reserve generation in GW", "techs":["CCGT", "CCGT_H2", "GT", "GT_H2"], "commodity":"Elec", "scaling":1, "cmap": "viridis"},
            }



                    
    cb = res["commodity_balance"]

    for i, j in map_dict.items():
        # Check if the data is being accessed correctly
        print("Accessing data from R?EMix model results...")
        try:
            data = cb.loc[idx[:,:,:,:,j["techs"],j["commodity"]], idx[:]]
        except Exception as e:
            print(f"Error accessing data: {e}")

        data = data.groupby("timeModel").sum().mul(j["scaling"]).values.reshape(-1, 24).T


        # Check if the figure is being created and saved correctly
        print("Creating and saving figure...")
        try:
            fig = plt.figure(figsize=(12.5,4), layout="constrained")
            ax1 = fig.add_subplot(111)
            # ... continue with plotting code ...

            # Check if the directory for saving figures exists
            Path(filename).parent.mkdir(parents=True, exist_ok=True)

            plt.savefig(f"{filename}.png" , dpi=300)
            plt.savefig(f"{filename}.svg")

        except Exception as e:
            print(f"Error creating or saving figure: {e}")



        fig = plt.figure(figsize=(12.5,4), layout="constrained")
        ax1 = fig.add_subplot(111)

        plt.imshow(data, aspect='auto', origin='upper', cmap=j["cmap"], interpolation='nearest', vmin=0,  vmax=np.max(data))
        cbar = plt.colorbar(pad=.02)
        if "text" in j.keys():
            ax1.text(-0.09, 0.96, j["text"], fontsize=16, transform=ax1.transAxes)

        plt.grid(False)

        # setting ticks positions
        t = np.arange(-0.5, 364.6, 30)
        ax1.xaxis.set_ticks(t)
        ax1.set_xticklabels(((t + 0.5)).astype(int))

        t = np.arange(-0.5, 23.6, 6)
        ax1.yaxis.set_ticks(t)
        ax1.set_yticklabels(list(map(lambda x: x + ":00", (t + 0.5).astype(int).astype(str))))

        cbar.set_label(j["title"], labelpad=15) # configure color bar

        filename = "{path_result}/figures/isopleths/iso_charging_{}".format(i)
        Path(filename).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(f"{filename}.png" , dpi=300)
        plt.savefig(f"{filename}.svg")

        plt.show()

carnot_isopleths_gen(res)


# %%