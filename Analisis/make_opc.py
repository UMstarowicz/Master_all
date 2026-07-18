from openmm.app import Modeller, Topology, PDBFile, ForceField
from openmm import unit, Vec3

topology = Topology()
positions = []

modeller = Modeller(topology, positions)

# Use TIP3P ONLY to generate geometry
forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")

modeller.addSolvent(
    forcefield,
    model="tip3p",
    boxSize=Vec3(3.0, 3.0, 3.0) * unit.nanometer
)

with open("water_opc.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

print("✅ water_opc.pdb created (geometry)")