from openmm import *
from openmm.app import *
from openmm.unit import *

# =========================
# PARAMETERS
# =========================
temperature = 295 * kelvin

timestep = 2 * femtoseconds
friction = 1 / picosecond

nsteps_equil = 500_000     # 1 ns
nsteps_prod  = 5_000_000   # 10 ns
report_interval = 5000

pdb_file = "water_tip4p2005.pdb"
output_prefix = "opc_295K_NVT"

# =========================
# LOAD STRUCTURE
# =========================
pdb = PDBFile(pdb_file)

# =========================
# FORCE FIELD (OPC – CORRECT)
# =========================
forcefield = ForceField(
    "opc.xml"
)

system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0 * nanometer,
    constraints=HBonds,
    rigidWater=True
)

# =========================
# INTEGRATOR (NVT)
# =========================
integrator = LangevinIntegrator(
    temperature,
    friction,
    timestep
)

# =========================
# PLATFORM (OpenCL)
# =========================
platform = Platform.getPlatformByName("OpenCL")
properties = {"OpenCLPrecision": "mixed"}

simulation = Simulation(
    pdb.topology,
    system,
    integrator,
    platform,
    properties
)

simulation.context.setPositions(pdb.positions)

# =========================
# MINIMIZATION
# =========================
print("Minimizing energy...")
simulation.minimizeEnergy()

# =========================
# EQUILIBRATION
# =========================
print("Equilibrating (NVT)...")
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(nsteps_equil)

# =========================
# REPORTERS
# =========================
simulation.reporters.append(
    DCDReporter(f"{output_prefix}.dcd", report_interval)
)

simulation.reporters.append(
    StateDataReporter(
        f"{output_prefix}.log",
        report_interval,
        step=True,
        time=True,
        temperature=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True
    )
)

from openmm.app import PDBFile, ForceField, Simulation, DCDReporter, StateDataReporter
from openmm.unit import *
from openmm import Platform, LangevinIntegrator
import numpy as np

# =========================
# Funkcja do zapisu XYZ
# =========================
class XYZReporter:
    """Reporter zapisujący trajektorie w formacie XYZ."""
    def __init__(self, file, reportInterval, topology):
        self._file = open(file, 'w')
        self._interval = reportInterval
        self._topology = topology
        self._n_atoms = self._topology.getNumAtoms()
        self._step = 0

    def describeNextReport(self, simulation):
        steps = self._interval - self._step % self._interval
        return (steps, True, False, False, False, False)

    def report(self, simulation, state):
        positions = state.getPositions(asNumpy=True) / nanometer  # w nm -> XYZ w nm
        self._file.write(f"{self._n_atoms}\n")
        self._file.write(f"Step {simulation.currentStep}\n")
        for atom, pos in zip(self._topology.atoms(), positions):
            self._file.write(f"{atom.element.symbol} {pos[0]:.5f} {pos[1]:.5f} {pos[2]:.5f}\n")
        self._file.flush()
        self._step = simulation.currentStep

    def __del__(self):
        self._file.close()


# =========================
# Dodanie XYZReporter
# =========================
xyz_file = f"{output_prefix}.xyz"
simulation.reporters.append(
    XYZReporter(xyz_file, report_interval, pdb.topology)
)


# =========================
# PRODUCTION
# =========================
print("Running production...")
simulation.step(nsteps_prod)

print("Simulation finished.")