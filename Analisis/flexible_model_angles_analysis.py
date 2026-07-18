import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# =========================
# LOAD DATA
# =========================
fname = "HOH_angles_amoeba_260K_NVT.dat"
angles = np.loadtxt(fname, skiprows=1)

angles = angles[~np.isnan(angles)]

mean = np.mean(angles)
std = np.std(angles)

print("Mean angle:", mean)
print("Std angle:", std)

# =========================
# HISTOGRAM
# =========================
plt.figure(figsize=(6,5))

# histogram (density=True = probability density)
counts, bins, _ = plt.hist(
    angles,
    bins=80,
    density=True,
    alpha=0.6
)

# Gaussian fit
x = np.linspace(min(angles), max(angles), 500)
gauss = norm.pdf(x, mean, std)

plt.plot(x, gauss)

plt.xlabel("HOH angle [deg]")
plt.ylabel("Probability density")
plt.title("AMOEBA HOH Angle Distribution")
plt.tight_layout()
plt.savefig("amoeba_HOH_histogram.png")
plt.show()
