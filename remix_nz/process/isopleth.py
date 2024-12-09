import pandas as pd
import seaborn as sns
import matplotlib as plt

sns.set(color_codes=True)
sns.set_style("whitegrid")
idx = pd.IndexSlice


def carnot_isopleths_gen(res):
    map_dict = {
                # "Electrolyzer": {"title":"Electrolyzer output in GW", "techs":"Electrolyzer", "commodity":"H2", "scaling":1, "cmap": "viridis"},
    #             "Methaniser": {"title":"Methanizer output in GW", "techs":"Methaniser", "commodity":"CH4", "scaling":1, "cmap": "viridis"},
                "PV": {"title":"PV generation in GW", "techs":"Photovoltaic", "commodity":"Elec", "scaling":1, "cmap": "viridis", "text": "b)"},
                "Wind": {"title":"Wind generation in GW", "techs":["WindOnshore","WindOffshore"], "commodity":"Elec", "scaling":1, "cmap": "viridis", "text": "a)"},
                # "CSP": {"title":"CSP generation in GW", "techs":"CSP_Powerblock", "commodity":"Elec", "scaling":1, "cmap": "viridis"},
    #            "RES": {"title":"Renewable generation in GW", "techs":["Photovoltaic","WindOnshore","WindOffshore","CSP_Powerblock"], "commodity":"Elec", "scaling":1, "cmap": "viridis"},
    #            "CHP": {"title":"CHP generation in GW", "techs":['TH_ExCCGT_NGas_XL','DH_Engine_NGas_M','DH_Engine_NGas_M', 'DH_FuelCell_H2_M'], "commodity":"Elec", "scaling":1, "cmap": "viridis"},
    #            "GT": {"title":"Reserve generation in GW", "techs":["CCGT", "CCGT_H2", "GT", "GT_H2"], "commodity":"Elec", "scaling":1, "cmap": "viridis"},
            }

    cb = res["commodity_balance"]

    for i, j in map_dict.items():
        data = cb.loc[idx[:,:,:,"2050",j["techs"],j["commodity"]], idx[:]]

        data = data.groupby("timeModel").sum().mul(j["scaling"]).values.reshape(-1, 24).T

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

        filename = "figures/isopleths/iso_charging_{}".format(i)
        Path(filename).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(f"{filename}.png" , dpi=300)
        plt.savefig(f"{filename}.svg")

        plt.show()

carnot_isopleths_gen()