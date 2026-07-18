import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from scipy.stats import skew, kurtosis

# =========================
# FIND ALL ANGLE FILES
# =========================
files = glob.glob("HOH_angles_*.dat")

results = {}

plt.figure(figsize=(7,5))

for fname in files:

    # ---- extract model info from filename ----
    base = os.path.basename(fname)
    label = base.replace("HOH_angles_", "").replace(".dat","")

    # ---- load angles ----
    angles = np.loadtxt(fname, skiprows=1)

    # ---- remove possible NaNs ----
    angles = angles[~np.isnan(angles)]

    # ---- statistics ----
    mean = np.mean(angles)
    std = np.std(angles)
    min_val = np.min(angles)
    max_val = np.max(angles)
    sk = skew(angles)
    kurt = kurtosis(angles)

    results[label] = {
        "Mean": mean,
        "Std": std,
        "Min": min_val,
        "Max": max_val,
        "Skewness": sk,
        "Kurtosis": kurt
    }

    # ---- histogram ----
    plt.hist(angles, bins=100, density=True, alpha=0.4, label=label)

plt.xlabel("HOH angle [deg]")
plt.ylabel("Probability density")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig("HOH_angle_distributions.png")
plt.show()

# =========================
# PRINT SUMMARY
# =========================
print("\nHOH Angle Structural Analysis\n")

for model, vals in results.items():
    print(f"{model:30s} "
          f"Mean={vals['Mean']:.3f}  "
          f"Std={vals['Std']:.3f}  "
          f"Min={vals['Min']:.2f}  "
          f"Max={vals['Max']:.2f}  "
          f"Skew={vals['Skewness']:.3f}")
