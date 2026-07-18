from openmm import *
from openmm.app import *
from openmm.unit import *

# =========================
# PARAMETERS
# =========================
temperature = 260 * kelvin
pressure = 1.0 * bar

timestep = 0.5 * femtoseconds   # wymagane dla Drude
friction = 1 / picosecond
drude_friction = 20 / picosecond

nsteps_equil = 500_000
nsteps_prod  = 5_000_000

report_interval = 5000
barostat_frequency = 25

pdb_file = "water_tip3p.pdb"
output_prefix = "swm4ndp_260 K_NPT"

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
# MODELLER
# =========================
modeller = Modeller(pdb.topology, pdb.positions)

# ADD DRUDE + M-SITE
modeller.addExtraParticles(forcefield)

# =========================
# SYSTEM
# =========================
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0 * nanometer,
    constraints=AllBonds,
    rigidWater=True
)

# NPT BAROSTAT
system.addForce(
    MonteCarloBarostat(
        pressure,
        temperature,
        barostat_frequency
    )
)

# =========================
# DRUDE INTEGRATOR
# =========================
integrator = DrudeLangevinIntegrator(
    temperature,
    friction,
    temperature,       # Drude temperature
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
print("Equilibrating (NPT)...")
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
        totalEnergy=True,
        volume=True,
        density=True
    )
)

# =========================
# PRODUCTION
# =========================
print("Running production (NPT)...")
simulation.step(nsteps_prod)

print("SWM4-DP NPT simulation finished")
