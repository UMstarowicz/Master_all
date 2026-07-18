import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

# =========================
# INPUT
# =========================
topology = "water_tip3p.pdb"
trajectory = "tip3p_295K_NVT.dcd"

r_max = 1.0      # nm
dr = 0.01        # nm
bins = int(r_max / dr)

# =========================
# LOAD SYSTEM
# =========================
u = mda.Universe(topology, trajectory)

oxygen = u.select_atoms("name O")
hydrogen = u.select_atoms("name H1 H2")

hist_OH = np.zeros(bins)
hist_HH = np.zeros(bins)

# =========================
# RDF LOOP
# =========================
for ts in u.trajectory:
    box = ts.dimensions[:3]

    # --- gOH ---
    for O in oxygen:
        for H in hydrogen:
            if O.resid == H.resid:
                continue
            r = np.linalg.norm(O.position - H.position)
            if r < r_max * 10:   # nm → Å
                bin_idx = int(r / (dr * 10))
                hist_OH[bin_idx] += 1

    # --- gHH ---
    for i, H1 in enumerate(hydrogen):
        for H2 in hydrogen[i+1:]:
            if H1.resid == H2.resid:
                continue
            r = np.linalg.norm(H1.position - H2.position)
            if r < r_max * 10:
                bin_idx = int(r / (dr * 10))
                hist_HH[bin_idx] += 2  # symmetry

# =========================
# NORMALIZATION
# =========================
rho_H = len(hydrogen) / u.dimensions[0]**3
rho_O = len(oxygen) / u.dimensions[0]**3

r = np.linspace(0, r_max, bins)

shell_vol = 4/3 * np.pi * ((r+dr)**3 - r**3)

g_OH = hist_OH / (len(oxygen) * rho_H * shell_vol * len(u.trajectory))
g_HH = hist_HH / (len(hydrogen) * rho_H * shell_vol * len(u.trajectory))

# =========================
# PLOT
# =========================
plt.figure()
plt.plot(r, g_OH, label="gOH")
plt.plot(r, g_HH, label="gHH")
plt.xlabel("r [nm]")
plt.ylabel("g(r)")
plt.legend()
plt.tight_layout()
plt.savefig("gOH_gHH.png")
plt.show()
