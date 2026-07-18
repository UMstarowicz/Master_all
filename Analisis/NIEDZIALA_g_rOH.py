import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt
import numpy as np

# =========================
# INPUT
# =========================
u = mda.Universe("water_tip3p.pdb", "tip3p_295K_NVT.dcd")
output_file = "gOO_tip3p_260K_NVT.dat"

O = u.select_atoms("name O")
H = u.select_atoms("name H1 H2 H")

# =========================
# RDF O–H (FIXED!)
# =========================
rdf_oh = InterRDF(
    O,
    H,
    range=(0.0, 6.0),
    nbins=200,
    exclude_same_residue=True   # <<< KLUCZOWE
)

rdf_oh.run()

r = rdf_oh.results.bins
g = rdf_oh.results.rdf

# =========================
# SAVE
# =========================
np.savetxt(
    output_file,
    np.column_stack((r, g)),
    header="r [A]   g_OH(r)",
    comments=""
)

# Plot
plt.plot(rdf_oh.results.bins, rdf_oh.results.rdf)
plt.xlabel("r (Å)")
plt.ylabel("g_OO(r)")
plt.title("O–H RDF, TIP3P, 260 K, NVT")
plt.show()

print("RDF O–H zapisany poprawnie")
