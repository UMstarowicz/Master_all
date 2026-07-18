from openmm import *
from openmm.app import *
from openmm.unit import *

# =========================
# PARAMETERS
# =========================
temperature = 295 * kelvin
pressure = 1.0 * bar

timestep = 2 * femtoseconds
friction = 1 / picosecond

nsteps_equil = 500_000     # 1 ns
nsteps_prod  = 5_000_000   # 10 ns
report_interval = 5000

pdb_file = "water_tip3p.pdb"
output_prefix = "tip3p_295K_NPT"

# =========================
# LOAD STRUCTURE
# =========================
pdb = PDBFile(pdb_file)

# =========================
# FORCE FIELD
# =========================
forcefield = ForceField("tip3p.xml")

system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0 * nanometer,
    constraints=HBonds,
    rigidWater=True
)

# =========================
# BAROSTAT (NPT)
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
# PLATFORM
# =========================
platform = Platform.getPlatformByName("CPU")

simulation = Simulation(
    pdb.topology,
    system,
    integrator,
    platform
)

simulation.context.setPositions(pdb.positions)

# =========================
# ENERGY MINIMIZATION
# =========================
print("Minimizing energy...")
simulation.minimizeEnergy()

# =========================
# EQUILIBRATION (NPT)
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