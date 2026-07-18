import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# SETTINGS
# =========================
equil_time_ps = 200.0

# =========================
# FILES
# =========================
model_map = {
    "opc_260K_NPT.log": "OPC",
    "tip4p2005_260K_NPT.log": "TIP4P/2005",
    "tip3p_260K_NPT.log": "TIP3P",
    "tip5p_260K_NPT.log": "TIP5P",
}

# =========================
# KDE DISTRIBUTIONS
# =========================
plt.figure(figsize=(10, 6))

for filepath, label in model_map.items():

    print(f"Processing {filepath}")

    df = pd.read_csv(filepath, comment="#", sep=",", header=None)

    df.columns = [
        "Step",
        "Time (ps)",
        "Potential Energy (kJ/mole)",
        "Kinetic Energy (kJ/mole)",
        "Total Energy (kJ/mole)",
        "Temperature (K)",
        "Pressure (bar)",
        "Density (g/cm3)",
    ]

    # Remove equilibration
    df_prod = df[df["Time (ps)"] >= equil_time_ps]

    temp = df_prod["Temperature (K)"]

    mean_T = temp.mean()

    # KDE curve
    line = temp.plot(
        kind="kde",
        linewidth=2.5,
        label=f"{label} (μ={mean_T:.2f} K)"
    )

    # Get the color used for KDE
    color = line.get_lines()[-1].get_color()

    # Mean line in same color
    plt.axvline(
        mean_T,
        linestyle="--",
        linewidth=2,
        color=color,
        alpha=0.9
    )

plt.xlabel("Temperature (K)", fontsize=12)
plt.ylabel("Probability density", fontsize=12)
plt.title("Temperature distribution in NPT ensemble", fontsize=14)

plt.legend(frameon=False)

# Cleaner publication style
plt.grid(alpha=0.2)
plt.tight_layout()

plt.savefig(
    "all_models_temperature_distribution_KDE_260K_NVP.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# =========================
# DENSITY KDE DISTRIBUTIONS
# =========================
plt.figure(figsize=(10, 6))

for filepath, label in model_map.items():

    print(f"Processing density: {filepath}")

    df = pd.read_csv(filepath, comment="#", sep=",", header=None)

    df.columns = [
        "Step",
        "Time (ps)",
        "Potential Energy (kJ/mole)",
        "Kinetic Energy (kJ/mole)",
        "Total Energy (kJ/mole)",
        "Temperature (K)",
        "Pressure (bar)",
        "Density (g/cm3)",
    ]

    # Remove equilibration
    df_prod = df[df["Time (ps)"] >= equil_time_ps]

    density = df_prod["Density (g/cm3)"]

    mean_rho = density.mean()

    # KDE curve
    line = density.plot(
        kind="kde",
        linewidth=2.5,
        label=f"{label} (μ={mean_rho:.4f})"
    )

    # Match mean line colour
    color = line.get_lines()[-1].get_color()

    plt.axvline(
        mean_rho,
        linestyle="--",
        linewidth=2,
        color=color,
        alpha=0.9
    )

plt.xlabel("Density (g/cm³)", fontsize=12)
plt.ylabel("Probability density", fontsize=12)
plt.title("Density distribution in NPT ensemble", fontsize=14)

plt.legend(frameon=False)

plt.grid(alpha=0.2)
plt.tight_layout()

plt.savefig(
    "all_models_density_distribution_KDE_260K_NPT.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()