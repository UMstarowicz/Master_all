from openmm.app import *
from openmm import *
from openmm.unit import *

# =========================
# INPUT
# =========================
input_pdb  = "water_tip3p.pdb"
output_pdb = "water_swm4ndp.pdb"

# =========================
# LOAD PDB
# =========================
pdb = PDBFile(input_pdb)

# =========================
# FORCE FIELD (SWM4-DP)
# =========================
forcefield = ForceField(
    "amber14-all.xml",
    "swm4ndp.xml"
)

# =========================
# MODELLER
# =========================
modeller = Modeller(pdb.topology, pdb.positions)

# KLUCZOWA LINIA
modeller.addExtraParticles(forcefield)

# =========================
# WRITE PDB
# =========================
with open(output_pdb, "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

print("SWM4-DP topology written:", output_pdb)
print("Number of atoms:", modeller.topology.getNumAtoms())
