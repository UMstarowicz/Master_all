from openmm import *
from openmm.app import *
from openmm.unit import *

# =========================
# PARAMETERS
# =========================
temperature = 260 * kelvin
timestep = 2 * femtoseconds
friction = 1 / picosecond
nsteps_equil = 500_000     # 1 ns equilibration
nsteps_prod = 5_000_000    # 10 ns production
report_interval = 5000

pdb_file = "water_tip3p.pdb"
output_prefix = "tip3p_260K_NVT"

# =========================
# LOAD STRUCTURE
# =========================
pdb = PDBFile(pdb_file)

# =========================
# FORCE FIELD (TIP3P)
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
# INTEGRATOR (NVT)
# =========================
integrator = LangevinIntegrator(
    temperature,
    friction,
    timestep
)

# =========================
# PLATFORM
# =========================
platform = Platform.getPlatformByName("CPU")  # change to "CPU" if needed

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
# EQUILIBRATION (NVT)
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
# PRODUCTION RUN
# =========================
print("Running production...")
simulation.step(nsteps_prod)

print("Simulation finished.")