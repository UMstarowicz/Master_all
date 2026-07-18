from openmm.app import Modeller, Topology, PDBFile, ForceField
from openmm import unit, Vec3

# 1. Pusta topologia
topology = Topology()
positions = []

modeller = Modeller(topology, positions)

# 2. Forcefield
forcefield = ForceField(
    "charmm36.xml",
    "charmm36/tip4p2005.xml"
)

# 3. Dodaj wodę
modeller.addSolvent(
    forcefield,
    model="tip4pew",
    boxSize=Vec3(3, 3, 3) * unit.nanometer
)

# 4. Zapis do PDB
with open("water_tip4p2005.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

print("✅ water_tip4p2005.pdb created")