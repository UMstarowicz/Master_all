import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# =========================
# FIND AMOEBA FILES
# =========================
files = sorted(glob.glob("HOH_angles_amoeba_*K_*.dat"))

# =========================
# LOAD ALL DATA FIRST
# =========================
all_angles = []
data_dict = {}

for fname in files:
    label = os.path.basename(fname).replace("HOH_angles_", "").replace(".dat","")
    angles = np.loadtxt(fname, skiprows=1)
    angles = angles[~np.isnan(angles)]

    data_dict[label] = angles
    all_angles.extend(angles)

all_angles = np.array(all_angles)

# =========================
# DEFINE COMMON BINS
# =========================
nbins = 120
bins = np.linspace(np.min(all_angles),
                   np.max(all_angles),
                   nbins)

# =========================
# PLOT
# =========================
plt.figure(figsize=(7,5))

for label, angles in data_dict.items():

    hist, edges = np.histogram(
        angles,
        bins=bins,
        density=True
    )

    centers = 0.5 * (edges[1:] + edges[:-1])

    plt.plot(centers, hist, label=label)

plt.xlabel("HOH angle [deg]")
plt.ylabel("Probability density")
plt.title("AMOEBA HOH Angle Distributions")
plt.legend()
plt.tight_layout()
plt.savefig("amoeba_HOH_comparison.png")
plt.show()
