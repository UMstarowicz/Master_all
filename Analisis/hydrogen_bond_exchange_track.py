import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

# =========================
# INPUT
# =========================
topology = "water_tip3p.pdb"
trajectory = "tip3p_295K_NVT.dcd"

# Hydrogen bond criteria
distance_cutoff = 3.5  # Å (O-O distance)
angle_cutoff = 150.0   # degrees

# =========================
# LOAD SYSTEM
# =========================
u = mda.Universe(topology, trajectory)

# =========================
# HYDROGEN BOND ANALYSIS
# =========================
hbond_analysis = HydrogenBondAnalysis(
    universe=u,
    donors_sel="name O",
    hydrogens_sel="name H1 H2",
    acceptors_sel="name O",
    d_a_cutoff=distance_cutoff,
    d_h_a_angle_cutoff=angle_cutoff
)

hbond_analysis.run()

# Results array columns:
# frame, donor index, hydrogen index, acceptor index, distance, angle
hbonds = hbond_analysis.results.hbonds

# =========================
# TRACK EXCHANGES
# =========================
frames = np.unique(hbonds[:, 0]).astype(int)
exchange_counts = []
previous_pairs = set()

for frame in frames:
    frame_hbonds = hbonds[hbonds[:, 0] == frame]

    # store donor-acceptor pairs (ignore hydrogen index)
    current_pairs = set((int(row[1]), int(row[3])) for row in frame_hbonds)

    # count changes compared to previous frame
    broken = previous_pairs - current_pairs
    formed = current_pairs - previous_pairs

    exchange_events = len(broken) + len(formed)
    exchange_counts.append(exchange_events)

    previous_pairs = current_pairs

# =========================
# SAVE DATA
# =========================
time_ps = frames * u.trajectory.dt  # convert frame index to time

data = np.column_stack((time_ps, exchange_counts))
np.savetxt("hbond_exchange.dat",
           data,
           header="Time(ps)    ExchangeEvents",
           fmt="%.6f")

# =========================
# PLOT
# =========================
plt.figure()
plt.plot(time_ps, exchange_counts)
plt.xlabel("Time [ps]")
plt.ylabel("H-bond Exchange Events")
plt.title("Hydrogen Bond Network Rearrangement")
plt.tight_layout()
plt.savefig("hbond_exchange.png")
plt.show()

# =========================
# SUMMARY METRICS
# =========================
avg_exchange_rate = np.mean(exchange_counts) / u.trajectory.dt

print("Average exchange events per frame:", np.mean(exchange_counts))
print("Approximate exchange rate (events/ps):", avg_exchange_rate)
