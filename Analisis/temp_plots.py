import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

log_file = "opc_295K_NVT.log"
output_prefix = "opc_295K_NVT"
equil_time_ps = 200.0

df = pd.read_csv(log_file)
df.columns = [c.strip() for c in df.columns]

print("Columns:", df.columns)

if df["Time (ps)"].max() < equil_time_ps:
    raise ValueError("Equilibration cutoff too large!")

df_prod = df[df["Time (ps)"] >= equil_time_ps]

print(f"Total frames: {len(df)}")
print(f"Production frames: {len(df_prod)}")

# ===== Stats =====
stats = {
    "Temperature mean [K]": df_prod["Temperature (K)"].mean(),
    "Temperature std [K]": df_prod["Temperature (K)"].std(ddof=0),
    "Temperature min [K]": df_prod["Temperature (K)"].min(),
    "Temperature max [K]": df_prod["Temperature (K)"].max(),
    "Potential energy mean [kJ/mol]": df_prod["Potential Energy (kJ/mole)"].mean(),
}

with open(f"{output_prefix}_summary.txt", "w") as f:
    for k, v in stats.items():
        f.write(f"{k}: {v:.6f}\n")

print(stats)

# ===== Plot function =====
def plot(x, y, xlabel, ylabel, fname, mean=None):
    plt.figure()
    plt.plot(x, y)
    if mean is not None:
        plt.axhline(mean, linestyle='--', label=f"Mean = {mean:.2f}")
        plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()

# ===== Plots =====
plot(df_prod["Time (ps)"], df_prod["Temperature (K)"],
     "Time (ps)", "Temperature (K)",
     f"{output_prefix}_temperature.png",
     mean=stats["Temperature mean [K]"])

plot(df_prod["Time (ps)"], df_prod["Potential Energy (kJ/mole)"],
     "Time (ps)", "Potential Energy (kJ/mol)",
     f"{output_prefix}_potential_energy.png")

print("Done.")