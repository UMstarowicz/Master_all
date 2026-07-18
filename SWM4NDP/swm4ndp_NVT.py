from openmm import *
from openmm.app import *
from openmm.unit import *

# =========================
# PARAMETERS
# =========================
temperature = 295 * kelvin

timestep = 0.5 * femtoseconds   # 2 fs = niestabilne dla Drude

friction = 1 / picosecond
drude_friction = 20 / picosecond

nsteps_equil = 500_000
nsteps_prod  = 5_000_000

report_interval = 5000

pdb_file = "water_tip3p.pdb"
output_prefix = "swm4ndp_295K_NVT"

# =========================
# LOAD PDB
# =========================
pdb = PDBFile(pdb_file)

# =========================
# FORCE FIELD
# =========================
forcefield = ForceField(
    "amber14-all.xml",
    "swm4ndp.xml"
)

# =========================
# MODELLER (CRITICAL)
# =========================
modeller = Modeller(pdb.topology, pdb.positions)

# THIS ADDS DRUDE + M-SITE
modeller.addExtraParticles(forcefield)

# =========================
# SYSTEM (ONLY ONCE!)
# =========================
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0 * nanometer,
    constraints=AllBonds,
    rigidWater=True
)

# =========================
# DRUDE INTEGRATOR
# =========================
integrator = DrudeLangevinIntegrator(
    temperature,
    friction,
    temperature,
    drude_friction,
    timestep
)

integrator.setMaxDrudeDistance(0.02)  # nm

# =========================
# PLATFORM
# =========================
platform = Platform.getPlatformByName("OpenCL")
properties = {"OpenCLPrecision": "mixed"}

simulation = Simulation(
    modeller.topology,
    system,
    integrator,
    platform,
    properties
)

simulation.context.setPositions(modeller.positions)

# =========================
# MINIMIZATION
# =========================
print("Minimizing energy...")
simulation.minimizeEnergy()

# =========================
# EQUILIBRATION
# =========================
print("Equilibrating...")
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

print("SWM4-DP simulation finished")
