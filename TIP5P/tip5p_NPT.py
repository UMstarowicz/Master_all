from openmm import *
from openmm.app import *
from openmm.unit import *

# =========================
# PARAMETERS
# =========================
temperature = 295   * kelvin
pressure = 1.0 * bar

timestep = 2 * femtoseconds
friction = 1 / picosecond

nsteps_equil = 500_000     # 1 ns
nsteps_prod  = 5_000_000   # 10 ns
report_interval = 5000

pdb_file = "water_tip5p.pdb"
output_prefix = "tip5p_295K_NPT_lol"

# =========================
# LOAD STRUCTURE
# =========================
pdb = PDBFile(pdb_file)

# =========================
# FORCE FIELD (TIP4P/2005)
# =========================
forcefield = ForceField(
    "charmm36.xml",
    "charmm36/tip5p.xml"
)

system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0 * nanometer,
    constraints=HBonds,
    rigidWater=True
)

# =========================
# BAROSTAT
# =========================
system.addForce(
    MonteCarloBarostat(
        pressure,
        temperature,
        25
    )
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
# GPU PLATFORM
# =========================
platform = Platform.getPlatformByName("OpenCL")
properties = {
    "OpenCLPrecision": "mixed"
}

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
print("Equilibrating (NPT)...")
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
        volume=True,
        density=True,
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

print("Simulation finished.")