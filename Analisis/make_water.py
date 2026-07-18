from openmm.app import *
from openmm import *
from openmm.unit import *

box_size = 3.0 #* nanometer

topology = Topology()
positions = []

forcefield = ForceField('tip4pew.xml')

modeller = Modeller(topology, positions)
modeller.addSolvent(
    forcefield,
    model='tip4pew',
    boxSize=Vec3(box_size, box_size, box_size)
)

with open('water_box.pdb', 'w') as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

print("Utworzono water_box.pdb")
