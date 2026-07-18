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

oxygen = u.select_atoms("name O")
hydrogen = u.select_atoms("name H1 H2")

# Exclude same-residue pairs
rdf = InterRDF(
    oxygen,
    hydrogen,
    nbins=nbins,
    range=(0.0, r_max),
    exclusion_block=(1, 2)   # exclude O-H within same residue
)

rdf.run()

r = rdf.results.bins
g_r = rdf.results.rdf

# =========================
# SAVE
# =========================
data = np.column_stack((r, g_r))
np.savetxt(f"gOH_{output}.dat", data,
           header="r (Angstrom)    gOH_inter",
           fmt="%.6f")

# =========================
# PLOT
# =========================
plt.figure()
plt.plot(r, g_r)
plt.xlabel("r [Å]")
plt.ylabel("gOH(r)")
plt.title(f"Intermolecular O-H RDF {output}")
plt.tight_layout()
plt.savefig(f"gOH_{output}.png")
plt.show()
