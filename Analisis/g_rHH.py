import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.rdf import InterRDF

# =========================
# INPUT
# =========================
output = "amoeba_295K_NVT"
topology = "water_tip3p.pdb"
trajectory = f"{output}.dcd"

r_max = 10.0   # Å
nbins = 200

# =========================
# LOAD SYSTEM
# =========================
u = mda.Universe(topology, trajectory)

hydrogen = u.select_atoms("name H1 H2")

# =========================
# RDF CALCULATION (INTERMOLECULAR ONLY)
# =========================
rdf = InterRDF(
    hydrogen,
    hydrogen,
    nbins=nbins,
    range=(0.0, r_max),
    exclude_same='residue'   # <-- THIS is the key fix
)

rdf.run()

r = rdf.results.bins
g_r = rdf.results.rdf

# =========================
# SAVE DATA
# =========================
data = np.column_stack((r, g_r))

np.savetxt(f"gHH_{output}.dat",
           data,
           header="r (Angstrom)    gHH_inter",
           fmt="%.6f")

# =========================
# PLOT
# =========================
plt.figure()
plt.plot(r, g_r)
plt.xlabel("r [Å]")
plt.ylabel("gHH(r)")
plt.title(f"Intermolecular H-H RDF {output}")
plt.tight_layout()
plt.savefig(f"gHH_{output}.png")
plt.show()
