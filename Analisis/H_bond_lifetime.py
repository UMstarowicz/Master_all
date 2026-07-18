import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

# =========================
# INPUT
# =========================
output = "tip5p_295K_NVT"
topology = "water_tip5p.pdb"
trajectory = f"{output}.dcd"

distance_cutoff = 3.5
angle_cutoff = 150.0

# =========================
# LOAD SYSTEM
# =========================
u = mda.Universe(topology, trajectory)

# =========================
# DETECT HYDROGEN BONDS
# =========================
h = HydrogenBondAnalysis(
    universe=u,
    donors_sel="name O",
    hydrogens_sel="name H1 H2",
    acceptors_sel="name O",
    d_a_cutoff=distance_cutoff,
    d_h_a_angle_cutoff=angle_cutoff
)

h.run()

# =========================
# COMPUTE LIFETIME
# =========================
tau, C_t = h.lifetime(intermittency=0)

# Convert frame index to time (ps)
time = tau * u.trajectory.dt

# =========================
# SAVE DATA
# =========================
data = np.column_stack((time, C_t))
np.savetxt(f"hbond_lifetime_{output}.dat",
           data,
           header="Time(ps)    C(t)",
           fmt="%.6f")

# =========================
# PLOT
# =========================
plt.figure()
plt.plot(time, C_t)
plt.xlabel("Time [ps]")
plt.ylabel("C(t)")
plt.title(f"Hydrogen Bond Autocorrelation {output}")
plt.tight_layout()
plt.savefig(f"hbond_lifetime_{output}_continuous.png")
plt.show()

# =========================
# INTEGRATED LIFETIME
# =========================
lifetime_ps = np.trapezoid(C_t, time)
print(f"Estimated hydrogen bond lifetime (ps) {output}:", lifetime_ps)
print(u.trajectory.dt)
print(len(time))
