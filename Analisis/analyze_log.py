import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# INPUT
# =========================
log_file = "opc_295K_NVT.log"
output_prefix = "opc_295K_NVT"
equil_time_ps = 200.0

# =========================
# LOAD DATA
# =========================
df = pd.read_csv(log_file)
df.columns = [c.strip() for c in df.columns]

print("Columns found:")
print(df.columns)

# =========================
# CUT EQUILIBRATION
# =========================
df_prod = df[df["Time (ps)"] >= equil_time_ps]

print(f"Loaded {len(df)} frames")
print(f"Using {len(df_prod)} frames after equilibration")

# =========================
# BASIC STATISTICS
# =========================
stats = {
    "Temperature mean [K]": df_prod["Temperature (K)"].mean(),
    "Temperature std [K]": df_prod["Temperature (K)"].std(),
    "Potential energy mean [kJ/mol]": df_prod["Potential Energy (kJ/mole)"].mean(),
}

if "Density (g/mL)" in df_prod.columns:
    stats["Density mean [g/mL]"] = df_prod["Density (g/mL)"].mean()
    stats["Density std [g/mL]"] = df_prod["Density (g/mL)"].std()
else:
    print("Density not present (NVT simulation)")


# =========================
# SAVE SUMMARY
# =========================
with open(f"{output_prefix}_summary.txt", "w") as f:
    for k, v in stats.items():
        f.write(f"{k}: {v:.6f}\n")

print("Summary written")

# =========================
# PLOTS
# =========================
def plot(x, y, xlabel, ylabel, fname):
    plt.figure()
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()

plot(df_prod["Time (ps)"], df_prod["Temperature (K)"],
     "Time (ps)", "Temperature (K)",
     f"{output_prefix}_temperature.png")

#plot(df_prod["Time (ps)"], df_prod["Density (g/mL)"],
#     "Time (ps)", "Density (g/mL)",
#     f"{output_prefix}_density.png")

plot(df_prod["Time (ps)"], df_prod["Potential Energy (kJ/mole)"],
     "Time (ps)", "Potential Energy (kJ/mol)",
     f"{output_prefix}_potential_energy.png")

print("Plots saved")
