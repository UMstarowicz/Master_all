from openmm import *
from openmm.app import *
from openmm.unit import *

# =========================
# PARAMETERS
# =========================
temperature = 277 * kelvin
timestep = 1.0 * femtoseconds   # AMOEBA: 1 fs
friction = 1 / picosecond

nsteps_equil = 500_000
nsteps_prod  = 5_000_000

report_interval = 5000

pdb_file = "water_tip3p.pdb"
output_prefix = "amoeba_277K_NVT"

# =========================
# LOAD PDB
# =========================
pdb = PDBFile(pdb_file)

# =========================
# AMOEBA FORCE FIELD
# =========================
forcefield = ForceField("amoeba2018.xml")

system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=0.9 * nanometer,
    constraints=AllBonds,
    polarization="mutual"
)

# =========================
# INTEGRATOR
# =========================
integrator = LangevinIntegrator(
    temperature,
    friction,
    timestep
)

# =========================
# PLATFORM
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

# =========================
# PRODUCTION
# =========================
print("Running production...")
simulation.step(nsteps_prod)

print("AMOEBA NVT finished.")
