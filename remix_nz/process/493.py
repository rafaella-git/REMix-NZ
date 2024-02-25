#%% [markdown]
# # Evaluate results for 493


# %%[markdown]
# ### Import dependencies
from remix.framework.tools.gdx import GDXEval
from remix.framework.tools.plots import plot_network, plot_choropleth
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from PIL import Image, ImageDraw
import os
from pathlib import Path
import warnings
# Suppress the FutureWarning
warnings.simplefilter(action='ignore', category=FutureWarning)
idx = pd.IndexSlice #often used shortcut

# %%[markdown]
# ### Define relevant lists to locate files
# Demand files available (different scenarios)
indx=3
files_lst=["medpop_evs_base","low_pop_out_base","med_pop_out_base","high_pop_out_base"] #493 as in the course 493, this is for Liv and Sam
files_lst=["medpop_evs","low_pop_out","med_pop_out","high_pop_out"] 

# Years run
yrs=[2020,2030,2040,2050] 
yrs=[2030,2040,2050]
yrs_str='-'.join([str(item) for item in yrs])
nodes_lst=["NIS","AKL","WTO","TRN","BOP","HBY","CEN","WEL","NEL","CAN","OTG"]

# %%[markdown]
# ### Define directories
path_base = "C:/Local/REMix"  
path_input = f"{path_base}/remix_nz/input"
path_output = f"{path_base}/remix_nz/output" 
path_geo = f"{path_input}/shapefiles"      
geofile="11regionsNZ.geojson"
shp_file = f"{path_geo}/11regionsNZ"
shp_attrcol = "id"


# case specific
demand_file=files_lst[indx]
case_name=f"{demand_file}_{yrs_str}"
path_result = f"{path_output}/{case_name}/result" 
figures_dir = Path(f"{path_output}/figures").mkdir(parents=True, exist_ok=True)
mini_fig_dir = Path(f"{path_result}/mini_figures").mkdir(parents=True, exist_ok=True)
path_figures = f"{path_output}/figures"
path_mini_fig =f"{path_result}/mini_figures"
gdx_file = f"{path_result}/{case_name}.gdx"
csv_file = f"{path_result}/{case_name}"
fig_file=f"{path_figures}/{case_name}"
minifig_file=f"{path_mini_fig}/{demand_file}"


# path_mini_fig #(folder)
# fig_file=f"{path_mini_fig}/{demand_file}"
# path_figures #(folder for collage only)
# for collage =f"{path_figures}/{case_name}.png"


# Figure set up
# define borders of the plot
lat = [-47.758239, -34.111702]
lon = [165.692103, 179.050919]
figsize_dual = (13.0, 6)
figsize_single = (7.5, 6)
figsize_dual_inv = (6, 11.0)
figsize_dual_inv_short = (4.2, 11.0)
figsize_single_inv = (6, 7.5)

technology_colors = {
    "CCGT": "#ff5f2d",
    "PV": "#ffc000",
    "Wind": "#9dc3e6",
    "HVDC": "#70ad47",
}
background_color = "#FFFFFF"
def update_fig_params():
    # params
    plt.rcParams.update({"figure.autolayout": True})  # use tight_layout
    plt.rcParams.update({"figure.titlesize": 20})  # size subtitle
    plt.rcParams.update({"figure.dpi": 75})  # size subtitle
    plt.rcParams.update({"savefig.dpi": 300})  # dpi
    plt.rcParams.update({"font.size": 16})  # default font size
    plt.rcParams.update({"axes.titlesize": 20})  # title size per subplot
    plt.rcParams.update({"axes.labelsize": 18})  # label size colormap
update_fig_params() 

# Recover results of the model runs
def get_caps(gams_file=gdx_file):
    if os.path.exists(gams_file):
        results = GDXEval(gams_file)  
        print(f"----Results file {case_name} exists!")  
        caps_gdx = results["converter_caps"]
        caps_gdx.to_csv(f"{csv_file}_caps.csv")
        #caps_gdx = caps_gdx[caps_gdx['value'] > 0.01].dropna() # remove all capacities with less than 10 MW
        # Replace 'techs' values based on mapping
        tech_list=list(set(caps_gdx.index.get_level_values('techs')))
        tech_map = {tech: "PV" if tech.startswith("pv") else "Wind" if tech.startswith("wind") else tech for tech in tech_list}
        # Reset the index
        caps_gdx = caps_gdx.reset_index()
        # Apply the mapping to the 'techs' column
        caps_gdx['techs'] = caps_gdx['techs'].map(tech_map)
        # Set the index back
        caps_gdx = caps_gdx.set_index(['accNodesModel', 'accYears', 'techs', 'commodity', 'capType'])
        caps_gdx = caps_gdx.groupby(caps_gdx.index.names).agg({'value': 'sum'})
        # Reset the index to have a flat structure
        caps_gdx.reset_index(inplace=True)
        caps_gdx.set_index(['accNodesModel', 'accYears', 'techs', 'commodity', 'capType'], inplace=True)
        index_labels = caps_gdx.index
        caps_gdx = caps_gdx.loc[idx[:, :, :, "Elec", "total"], :].round(2)
        caps_gdx.to_csv(f"{csv_file}_caps_grouped.csv")
        return caps_gdx
        # Until here the dataframe works
        # -------------------------------------
    else:
        print("Oops: the result file {case_name} does not exist.")
        return None

# (Old because ugly) Plot installed capacities in Choropleth Map  
def old_plot_caps(df_caps=caps,year='2050',case=demand_file):
    df_caps = df_caps.loc[idx[nodes_lst, f"{year}", :, "Elec", "total"], :].round(2)
    # separate info per tech
    map_pv = df_caps.loc[idx[:, :, "PV"], idx[:]].groupby("accNodesModel").sum()#.div(1e3)
    map_wind = df_caps.loc[idx[:, :, "Wind"], idx[:]].groupby("accNodesModel").sum()#.div(1e3)
    map_battery = df_caps.loc[idx[:, :, "Battery"], idx[:]].groupby("accNodesModel").sum()#.div(1e3)
    # define borders of the plot
    lat = [-47.758239, -34.111702]
    lon = [165.692103, 179.050919]
    # params
    plt.rcParams.update({"figure.autolayout": True})  # use tight_layout
    plt.rcParams.update({"figure.titlesize": 20})  # size subtitle
    plt.rcParams.update({"figure.dpi": 75})  # size subtitle
    plt.rcParams.update({"savefig.dpi": 300})  # dpi
    plt.rcParams.update({"font.size": 16})  # default font size
    plt.rcParams.update({"axes.titlesize": 20})  # title size per subplot
    plt.rcParams.update({"axes.labelsize": 18})  # label size colormap
    figsize_dual = (13.0, 6)
    figsize_single = (7.5, 6)
    # plot
    fig = plt.figure(figsize=figsize_dual)
    fig.suptitle(f"{year}", fontsize=22, ha='center')
    fig.text(0.5, 0.01, f"Demand file: {case}", fontsize=10, ha='center')
    fig.patch.set_facecolor(background_color)
    ax1 = fig.add_subplot(121)
    ax1.set_facecolor("#F0F8FF")
    plot_choropleth(
        map_pv,   #df
        shp_file,
        shp_attrcol,
        lat,
        lon,
        title="Installed PV capacity",
        clabel="Power in GW",
        cmap="Oranges",
        ax=ax1,
    )

    ax2 = fig.add_subplot(122)
    ax2.set_facecolor("#F0F8FF")
    plot_choropleth(
        map_wind,
        shp_file,
        shp_attrcol,
        lat,
        lon,
        title="Installed wind capacity",
        clabel="Power in GW",
        cmap="Blues",
        ax=ax2,
    )

# Group images in collage
def old_group_imgs(fig_name=minifig_file,w=1,h=3):
    # Create a list of image filenames
    image_files = [f"{fig_name}_cap_2030.png", f"{fig_name}_cap_2040.png", f"{fig_name}_cap_2050.png"]

    # Set the dimensions of the collage
    collage_width = 13*300*w  # Width of the collage in pixels
    collage_height = 6*300*h  # Height of the collage in pixels

    # Create a new blank image with the specified dimensions
    collage = Image.new("RGB", (collage_width, collage_height), (255, 255, 255))
    # Create a drawing context to paste images onto the collage
    draw = ImageDraw.Draw(collage)
    # Set the position to start pasting images
    x = 0
    y = 0
    # Paste each image into the collage
    for image_file in image_files:
        image = Image.open(image_file)
        collage.paste(image, (x, y))
        y += image.height  # Move down to the next row
    # Save the collage as a new image file
    collage.save(f"{fig_file}_cap.png")

# Finds the maximum value in the caps df
def max_cap(tech,df=caps,nodes=nodes_lst):
    sliced_df=df.loc[idx[nodes_lst, :, f"{tech}", "Elec", "total"],  idx[:]].round(2)#.div(1e3)
    max=round(sliced_df['value'].max(),2)
    return(max)



#%%
def plot_caps_h(df_caps,year,case=demand_file):

    # slice capacities dataframe for specified techs
    map_pv=df_caps.loc[idx[nodes_lst, f"{year}",  "PV", "Elec", "total"],  idx[:]].groupby("accNodesModel").sum().round(2)
    map_wind=df_caps.loc[idx[nodes_lst, f"{year}", "Wind", "Elec", "total"],  idx[:]].groupby("accNodesModel").sum().round(2)
    #map_pv.to_csv(f"{csv_file}_dfview_pv.csv")
    #map_wind.to_csv(f"{csv_file}_dfview_wind.csv")

    # overall figure
    fig = plt.figure(figsize=figsize_dual)
    fig.patch.set_facecolor(background_color)
    #fig.suptitle("Installed capacity", fontsize=24, ha='center')

    # subfigure ax1
    ax1 = fig.add_subplot(121)
    #ax1.set_facecolor("#F0F8FF")
    plot_choropleth(
        map_pv,   #df
        shp_file,
        shp_attrcol,
        lat,
        lon,
        #title=f"PV {year}",
        clabel="Power in GW",
        cmap="Oranges",
        ax=ax1, 
        limits=[0,max_pv*1.01]
        #limits=[0,4.5*1.001] 
    )
    ax1.text(0.02, 0.95, f"  {year}", fontsize=22, ha='left', va='top', color='gray', transform=ax1.transAxes)
    ax1.text(0.95, 0.01, f"{case}", fontsize=12, ha='right', va='bottom', color='gray', transform=ax1.transAxes)


    ax2 = fig.add_subplot(122)
    #ax2.set_facecolor("#F0F8FF")
    plot_choropleth(
        map_wind,
        shp_file,
        shp_attrcol,
        lat,
        lon,
        #title=f"Wind {year}",
        clabel="Power in GW",
        cmap="Blues",
        ax=ax2,
        limits=[0,max_wind*1.03] #1.03
        #limits=[0,0.45*1.1] 

    )
    ax2.text(0.02, 0.95, f"  {year}", fontsize=22, ha='left', va='top', color='gray', transform=ax2.transAxes)
    ax2.text(0.95, 0.01, f"{case}", fontsize=12, ha='right', va='bottom', color='gray', transform=ax2.transAxes)
    # define each figure directoy
    mini_fig=f"{minifig_file}_cap_h_{year}.png"
    plt.savefig(mini_fig)


#%%
#hola
def plot_caps_v(df_caps,year,bar_bol=False, case=demand_file):
    # slice capacities dataframe for specified techs
    map_pv=df_caps.loc[idx[nodes_lst, f"{year}",  "PV", "Elec", "total"],  idx[:]].groupby("accNodesModel").sum().round(2)
    map_wind=df_caps.loc[idx[nodes_lst, f"{year}", "Wind", "Elec", "total"],  idx[:]].groupby("accNodesModel").sum().round(2)
    # overall figure
    fig = plt.figure(figsize=figsize_dual_inv_short)
    if bar_bol: #modify to make them match
        fig = plt.figure(figsize=figsize_dual_inv)
    fig.patch.set_facecolor(background_color)
    # subfigure ax1
    ax1 = fig.add_subplot(211)
    #ax1.set_facecolor("#F0F8FF")
    plot_choropleth(
        map_pv,   #df
        shp_file,
        shp_attrcol,
        lat,
        lon,
        #title=f"PV {year}",
        clabel="Photovoltaic capacity in GW",
        cmap="Oranges",
        cbar=bar_bol,
        ax=ax1, 
        limits=[0,4.48*1.01]
        #limits=[0,max_pv*1.01]
    )
    ax2 = fig.add_subplot(212)
    #ax2.set_facecolor("#F0F8FF")
    plot_choropleth(
        map_wind,
        shp_file,
        shp_attrcol,
        lat,
        lon,
        #title=f"Wind {year}",
        clabel="Wind capacity in GW",
        cmap="Blues",
        cbar=bar_bol,
        ax=ax2,
        limits=[0,0.35*1.03]        
        #limits=[0,max_wind*1.03] #1.03
    )
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    ax1.text(0.02, 0.95, f"  {year}", fontsize=22, ha='left', va='top', color='gray', transform=ax1.transAxes)
    ax1.text(0.95, 0.01, f"{case}", fontsize=12, ha='right', va='bottom', color='gray', transform=ax1.transAxes)
    ax2.text(0.02, 0.95, f"  {year}", fontsize=22, ha='left', va='top', color='gray', transform=ax2.transAxes)
    ax2.text(0.95, 0.01, f"{case}", fontsize=12, ha='right', va='bottom', color='gray', transform=ax2.transAxes)
    # define each figure directoy
    mini_fig=f"{minifig_file}_cap_v_{year}.png"
    plt.savefig(mini_fig, bbox_inches='tight', pad_inches=0.1)



#%%

# Group images in collage
def group_imgs(fig_name=minifig_file):
    # Create a list of image filenames
    image_files = [f"{fig_name}_cap_v_2030.png", f"{fig_name}_cap_v_2040.png", f"{fig_name}_cap_v_2050.png"]

    # Get the maximum height among the images
    max_height = max(Image.open(image_file).height for image_file in image_files)

    # Calculate the width of the collage by summing the individual image widths
    collage_width = sum(Image.open(image_file).width for image_file in image_files)

    # Create a new blank image with the specified dimensions
    collage = Image.new("RGB", (collage_width, max_height), (255, 255, 255))

    # Set the initial position to start pasting images
    x = 0

    # Paste each image into the collage side by side
    for image_file in image_files:
        image = Image.open(image_file)
        collage.paste(image, (x, 0))
        x += image.width  # Move to the right for the next image
    
    # Save the collage as a new image file
    collage.save(f"{fig_name}_cap.png")
    # it also saves it to the main figures folder in output>figures
    collage.save(f"{fig_file}_cap.png")



# Define the dataframe with the capacities and export to CSV
caps = get_caps()
# Get the limits of the plots to make them uniform
max_pv=max_cap("PV",caps)
max_wind=max_cap("Wind",caps)
# plot horizontal: pv left and wind right
plot_caps_h(caps,'2030')
plot_caps_h(caps,'2040')
plot_caps_h(caps,'2050')
# plot vertical: pv up and wind down
plot_caps_v(caps,'2030')
plot_caps_v(caps,'2040')
plot_caps_v(caps,'2050',True)
group_imgs()


#%%

def create_gif(files_lst=files_lst, name="cap", duration=600):
    """
    Create an animated GIF from a list of image files.

    Args:
        files_lst (list): List of demand files, used to define image file paths.
        output_filename (str): Desired filename for the resulting GIF.
        duration (int, optional): Duration (in milliseconds) between frames.
    """
    # Create a list to store the image frames
    frames = []
    
    # Load and append each image to the frames list
    for files in files_lst:
        image_file=f"{path_figures}/{files}_{yrs_str}_cap.png"
        img = Image.open(image_file)
        frames.append(img)

    # Set the path for the output GIF
    output_gif_path = f"{path_figures}/{name}.gif"

    # Save the frames as an animated GIF
    frames[0].save(output_gif_path, save_all=True, append_images=frames[1:], loop=0, duration=duration)

    print(f"GIF created and saved at: {output_gif_path}")
create_gif()

# %%
