import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# =========================
# PARAMETERS
# =========================
equil_time_ps = 200

models = ["TIP3P", "TIP4P2005", "TIP5P", "OPC"]#, "SWM4NDP"]
temps = ["260K", "277K", "295K"]
ensembles = ["NPT"]

# =========================
# HELPER FUNCTION
# =========================
def load_data(model, temp, ensemble):
    fname = f"{model}_{temp}_{ensemble}.log"
    if not os.path.exists(fname):
        print(f"Missing: {fname}")
        return None

    # ---- read header safely ----
    with open(fname, 'r', encoding='utf-8') as f:
        header_line = f.readline()

    # remove everything before first quote
    header_line = header_line[header_line.find('"'):]

    # split and clean column names
    columns = [c.strip().strip('"') for c in header_line.split(',')]

    # ---- read numeric data ----
    df = pd.read_csv(
        fname,
        skiprows=1,
        names=columns
    )

    # convert to numeric
    df = df.apply(pd.to_numeric, errors='coerce')

    # DEBUG (run once if needed)
    # print(columns)
    # print(df.columns)

    if "Time (ps)" not in df.columns:
        print(f"Problem with columns in {fname}")
        print(df.columns)
        return None

    df = df[df["Time (ps)"] >= equil_time_ps]

    return df
# =========================
# HELPER FUNCTION: running average
# =========================
window =100
def running_average(series, window=window):
    return series.rolling(window=window, min_periods=1, center=True).mean()

# =========================
# GLOBAL PLOT with smoothed data and average line
# =========================
plt.figure(figsize=(8,6))

for model in models:
    for temp in temps:
        for ensemble in ensembles:
            df = load_data(model, temp, ensemble)
            if df is None:
                continue

            label = f"{model}-{temp}-{ensemble}"

            if ensemble == "NPT" and "Density (g/mL)" in df.columns:
                y = df["Density (g/mL)"]
                y_smooth = running_average(y, window=window)  # smoothed density
                avg_density = y.mean()  # overall time average
                plt.plot(df["Time (ps)"], y_smooth, label=label, alpha=0.6)
                plt.axhline(avg_density, color='red', linestyle='--', alpha=0.5)
            else:
                y = df["Temperature (K)"]
                y_smooth = running_average(y, window=window)
                plt.plot(df["Time (ps)"], y_smooth, label=label, alpha=0.6)

plt.xlabel("Time (ps)")
plt.ylabel("Density (g/mL) or Temperature (K)")
plt.legend(fontsize=6)
plt.tight_layout()
plt.savefig("GLOBAL_all_models_smoothed_v2.png")
plt.close()

import matplotlib.pyplot as plt
from cycler import cycler

# =========================
# 6 PLOTS: temp × ensemble with smoothed data + mean line in same color
# =========================
for temp in temps:
    for ensemble in ensembles:
        plt.figure(figsize=(7,5))
        ax = plt.gca()  # pobieramy aktualną oś

        # Automatyczne kolory jak w legendzie
        color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

        for i, model in enumerate(models):
            df = load_data(model, temp, ensemble)
            if df is None:
                continue

            if ensemble == "NPT":
                y = df["Density (g/mL)"]
                ylabel = "Density (g/mL)"
            else:
                y = df["Temperature (K)"]
                ylabel = "Temperature (K)"

            y_smooth = running_average(y, window=window)

            color = color_cycle[i % len(color_cycle)]  # przypisanie koloru
            ax.plot(df["Time (ps)"], y_smooth, label=model, color=color)
            avg_value = y.mean()
            ax.axhline(avg_value, color=color, linestyle='--', alpha=0.5)  # średnia w tym samym kolorze

        plt.xlabel("Time (ps)")
        plt.ylabel(ylabel)
        plt.title(f"{temp} - {ensemble}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{temp}_{ensemble}_all_models_smoothed_v2.png")
        plt.close()


# =========================
# Ensemble comparison per temperature with smoothed data + mean line in same color
# =========================
for temp in temps:
    plt.figure(figsize=(7,5))
    ax = plt.gca()
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for i, model in enumerate(models):
        for j, ensemble in enumerate(ensembles):
            df = load_data(model, temp, ensemble)
            if df is None:
                continue

            if ensemble == "NPT":
                y = df["Density (g/mL)"]
            else:
                y = df["Temperature (K)"]

            y_smooth = running_average(y, window=window)
            color = color_cycle[i % len(color_cycle)]  # kolor zależy od modelu
            ax.plot(df["Time (ps)"], y_smooth, label=f"{model}-{ensemble}", color=color)
            avg_value = y.mean()
            ax.axhline(avg_value, color=color, linestyle='--', alpha=0.5)

    plt.xlabel("Time (ps)")
    plt.ylabel("Density (NPT)")
    plt.title(f"Ensemble comparison at {temp}")
    plt.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(f"{temp}_ensemble_comparison_smoothed_v3.png")
    plt.close()
