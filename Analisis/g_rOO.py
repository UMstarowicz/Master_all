import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt
import numpy as np

# Load system
u = mda.Universe("water_tip3p.pdb", "amoeba_260K_NVT.dcd")

# =========================
# SAVE TO FILE
# =========================
output_file = "gOO_amoeba_260K_NVT.dat"

# Select oxygen atoms
O = u.select_atoms("name O")

# RDF O-O
rdf = InterRDF(
    O,
    O,
    range=(0.0, 10.0),     # Å
    nbins=200,
    exclusion_block=(1,1) # CRITICAL for water
)

rdf.run()

r = rdf.results.bins
g_r = rdf.results.rdf

np.savetxt(
    output_file,
    np.column_stack((r, g_r)),
    header="r [A]    g_OO(r)",
    comments=""
)

# Plot
plt.plot(rdf.results.bins, rdf.results.rdf)
plt.xlabel("r (Å)")
plt.ylabel("g_OO(r)")
plt.title("O–O RDF, amoeba, 260 K, NVT")
plt.show()
