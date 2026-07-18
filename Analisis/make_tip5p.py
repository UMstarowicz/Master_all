from openmm.app import Modeller, Topology, PDBFile, ForceField
from openmm import unit, Vec3

# =========================
# EMPTY TOPOLOGY
# =========================
topology = Topology()
positions = []

modeller = Modeller(topology, positions)

# =========================
# FORCE FIELD (TIP5P)
# =========================
forcefield = ForceField(
    "charmm36.xml",
    "charmm36/tip5p.xml"
)

# =========================
# ADD WATER BOX
# =========================
modeller.addSolvent(
    forcefield,
    model="tip5p",
    boxSize=Vec3(3.0, 3.0, 3.0) * unit.nanometer
)

# =========================
# WRITE PDB
# =========================
with open("water_tip5p.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

print("✅ water_tip5p.pdb created")