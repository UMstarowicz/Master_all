import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# =========================
# INPUT FILES
# =========================
files = {
    "TIP3P": "gOO_tip3p_295K_NVT.dat",
    "TIP4P/2005": "gOO_tip4p2005_295K_NVT.dat",
    "TIP5P": "gOO_tip5p_295K_NVT.dat",
    "OPC": "gOO_opc_295K_NVT.dat"#,
    #"SWM4NDP": "gOO_swm4ndp_295K_NVT.dat",
    #"AMOEBA": "gOO_amoeba_295K_NVT.dat"
}

# ===== density for coordination number =====
density = 0.997  # g/cm^3
molar_mass = 18.01528  # g/mol
NA = 6.02214076e23

# convert to number density [Å^-3]
rho = density * NA / molar_mass / 1e24

results = {}

plt.figure()

for model, fname in files.items():
    data = np.loadtxt(fname, skiprows=1)
    r = data[:,0]
    g = data[:,1]

    plt.plot(r, g, label=model)

    # -------------------------
    # FIRST PEAK
    # -------------------------
    peaks, _ = find_peaks(g, height=1.5)

    if len(peaks) == 0:
        continue

    p = peaks[0]
    r1 = r[p]
    g1 = g[p]

    # -------------------------
    # FIRST MINIMUM
    # -------------------------
    after_peak = g[p+1:]
    min_index_local = np.argmin(after_peak)
    m = p + 1 + min_index_local

    r_min = r[m]
    g_min = g[m]

    # -------------------------
    # AREA under first peak
    # -------------------------
    area = np.trapezoid(g[:m], r[:m])

    # -------------------------
    # COORDINATION NUMBER
    # -------------------------
    integrand = 4*np.pi*r[:m]**2 * g[:m]
    coordination = rho * np.trapezoid(integrand, r[:m])

    results[model] = {
        "r1 [Å]": r1,
        "g1": g1,
        "r_min [Å]": r_min,
        "g_min": g_min,
        "Area_peak": area,
        "Coordination": coordination
    }

plt.xlabel("r [Å]")
plt.ylabel("gOO(r)")
plt.legend()
plt.tight_layout()
plt.savefig("gOO_comparison_v2.png")
plt.show()

# =========================
# PRINT SUMMARY
# =========================
print("\nStructural comparison:\n")

for model, vals in results.items():
    print(f"{model:15s} "
          f"r1={vals['r1 [Å]']:.3f} Å  "
          f"g1={vals['g1']:.3f}  "
          f"r_min={vals['r_min [Å]']:.3f} Å  "
          f"Area={vals['Area_peak']:.3f}  "
          f"CN={vals['Coordination']:.3f}")
