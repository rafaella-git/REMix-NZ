# %% [markdown]
# # 🔍 Peak-Normalized Load Profile Analysis
# This notebook loads electricity demand data, normalizes each year by its peak demand (100% = 1.0),
# and visualizes patterns over time with all values expressed as % of yearly peak.

# %% [markdown]
# ## 📦 Imports and Setup
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px

# %% [markdown]
# ## 📂 Load Data and Normalize by Peak
# Load normalized demand profile and convert to peak-normalized (100% = yearly max)
df_normalized = pd.read_csv(
    "C:/Local/REMix/remix_nz/process/data-for-paper/Normalized-Load-Profile.csv",
    index_col=0
)
df_normalized.columns = df_normalized.columns.astype(int)

# Create peak-normalized version (each year's max = 1.0)
df_peak_normalized = df_normalized.div(df_normalized.max(axis=0), axis=1)
print("Verification - Max values per year:")
print(df_peak_normalized.max().round(3))

# %% [markdown]
# ## ⚙️ Configuration
short_mode = False  # Toggle faster rendering
mode_name = "short" if short_mode else "full"
export_path = f"figures/{mode_name}"
os.makedirs(export_path, exist_ok=True)

years_to_plot = [2020, 2030, 2050, 2070, 2100]
df_plot = df_peak_normalized.iloc[::20] if short_mode else df_peak_normalized
x_ticks = np.linspace(0, len(df_plot.index) - 1, 10).astype(int)
x_labels = [df_plot.index[i] for i in x_ticks]

# %% [markdown]
# ## 🔥 Heatmaps (Hourly x Daily)
def reshape_to_daily(series):
    """Reshape hourly data to daily (365x24) format"""
    values = series.values.reshape((365, 24))
    return pd.DataFrame(values, 
                       index=[f"Day {i+1}" for i in range(365)],
                       columns=[f"{h:02d}:00" for h in range(24)])

fig, axs = plt.subplots(1, 2, figsize=(18, 6))
for i, year in enumerate([2020, 2050]):
    data = reshape_to_daily(df_peak_normalized[year])
    sns.heatmap(data, ax=axs[i], cmap='viridis')
    axs[i].set_title(f"Hourly Heatmap - {year} (Peak-Normalized)")
    axs[i].set_xlabel("Hour of Day")
    axs[i].set_ylabel("Day of Year")
plt.tight_layout()
plt.savefig(f"{export_path}/heatmaps_peak_normalized.png", dpi=300)
plt.show()

# %% [markdown]
# ## 📈 Peak-Normalized Hourly Profiles
plt.figure(figsize=(14, 6))
for year in years_to_plot:
    plt.plot(df_plot.index, df_plot[year], label=str(year), alpha=0.8)
plt.title("Peak-Normalized Hourly Demand Profiles")
plt.xlabel("Hour of the Year")
plt.ylabel("Demand (% of yearly peak)")
plt.xticks(ticks=x_ticks, labels=x_labels, rotation=45)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(f"{export_path}/hourly_profiles_peak_normalized.png", dpi=300)
plt.close()

# %% [markdown]
# ## ❄️ Summer vs Winter (Peak-Normalized)
summer = df_peak_normalized.loc["t0433":"t0496"]  # Summer week
winter = df_peak_normalized.loc["t4345":"t4408"]  # Winter week

plt.figure(figsize=(10, 5))
plt.plot(summer.index, summer[2020], label="Summer", linestyle='-')
plt.plot(winter.index, winter[2020], label="Winter", linestyle='--')
plt.title("Summer vs Winter Demand (2020, Peak-Normalized)")
plt.xlabel("Hour of the Week")
plt.ylabel("Demand (% of yearly peak)")
plt.xticks(ticks=np.linspace(0, 167, 10), 
           labels=[f"{int(x)}" for x in np.linspace(0, 167, 10)])
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(f"{export_path}/summer_vs_winter_peak_normalized.png", dpi=300)
plt.close()

# %% [markdown]
# ## 🔄 Delta: 2050 vs 2020 (Peak-Normalized)
delta = df_peak_normalized[2050] - df_peak_normalized[2020]
delta_plot = delta.iloc[::20] if short_mode else delta

plt.figure(figsize=(12, 5))
plt.plot(delta_plot.index, delta_plot.values, color='crimson')
plt.axhline(0, color='gray', linestyle='--')
plt.title("Change in Peak-Normalized Demand: 2050 vs 2020")
plt.xlabel("Hour of the Year")
plt.ylabel("Δ Demand (% points)")
plt.xticks(ticks=x_ticks, labels=x_labels, rotation=45)
plt.grid(True)
plt.tight_layout()
plt.savefig(f"{export_path}/delta_2050_vs_2020_peak_normalized.png", dpi=300)
plt.close()

# %% [markdown]
# ## 📉 Load Duration Curve (Peak-Normalized)
plt.figure(figsize=(10, 5))
for year in [2020, 2050]:
    ldc = df_peak_normalized[year].sort_values(ascending=False).reset_index(drop=True)
    ldc = ldc.iloc[::20] if short_mode else ldc
    plt.plot(ldc.index, ldc.values, label=str(year))
plt.title("Load Duration Curve (Peak-Normalized)")
plt.xlabel("Hour Rank")
plt.ylabel("Demand (% of yearly peak)")
plt.legend()
plt.tight_layout()
plt.savefig(f"{export_path}/load_duration_curve_peak_normalized.png", dpi=300)
plt.close()

# %% [markdown]
# ## 🌐 Interactive Plotly Visualization
fig = go.Figure()
for year in years_to_plot:
    fig.add_trace(go.Scatter(
        x=df_peak_normalized.index,
        y=df_peak_normalized[year],
        mode='lines',
        name=str(year)
    ))
fig.update_layout(
    title="Peak-Normalized Demand Profiles",
    xaxis_title="Hour of the Year",
    yaxis_title="Demand (% of yearly peak)",
    template="plotly_white",
    hovermode="x unified"
)
fig.write_html(f"{export_path}/interactive_peak_normalized.html")

# %% [markdown]
# ## 📊 First 48 Hours (Peak-Normalized)
df_first_48 = df_peak_normalized.iloc[:48].copy()
df_first_48.index = range(1, 49)

plt.figure(figsize=(10, 5))
for year in [2020, 2050]:  # Plot selected years for clarity
    plt.plot(df_first_48.index, df_first_48[year], label=str(year))
plt.title("First 48 Hours of Year (Peak-Normalized)")
plt.xlabel("Hour")
plt.ylabel("Demand (% of yearly peak)")
plt.legend()
plt.tight_layout()
plt.savefig(f"{export_path}/first_48h_peak_normalized.png", dpi=300)
plt.close()