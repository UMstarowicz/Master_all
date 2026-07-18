import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from scipy.spatial import cKDTree
print("Script started")
# =========================
# INPUT
# =========================
topology = "water_tip3p.pdb"
trajectory = "amoeba_260K_NPT.dcd"

output_prefix = "analysis"

# H-bond definition
d_a_cutoff = 3.5   # Å
angle_cutoff = 150 # degrees

# =========================
# LOAD SYSTEM
# =========================
u = mda.Universe(topology, trajectory)
print("Number of frames:", len(u.trajectory))
O = u.select_atoms("name O")
H = u.select_atoms("name H1 H2")
for ts in u.trajectory:
    print("Frame:", ts.frame)
    break
print("Number of frames:", len(u.trajectory))
# =========================
# H-BOND ANALYSIS
# =========================
h = HydrogenBondAnalysis(
    universe=u,
    donors_sel="name O",
    hydrogens_sel="name H1 H2",
    acceptors_sel="name O",
    d_a_cutoff=d_a_cutoff,
    d_h_a_angle_cutoff=angle_cutoff
)

h.run()

hbonds = h.results.hbonds
print("Number of frames:", len(u.trajectory))
# =========================
# STORAGE
# =========================
q_list = []
na_list = []
w4_list = []
qcos_list = []

# =========================
# MAIN LOOP
# =========================
for ts in u.trajectory:

    O_pos = O.positions

    # ---------- KDTree for neighbors ----------
    tree = cKDTree(O_pos)
    neighbors = tree.query_ball_tree(tree, r=3.5)

    # ---------- TETRAHEDRAL ORDER q ----------
    q_frame = []

    for i, neigh in enumerate(neighbors):

        if len(neigh) < 5:
            continue

        # take 4 nearest neighbors
        neigh = neigh[1:5]
        vecs = O_pos[neigh] - O_pos[i]

        cos_angles = []
        for j in range(3):
            for k in range(j+1, 4):
                v1 = vecs[j]
                v2 = vecs[k]
                cos = np.dot(v1, v2) / (np.linalg.norm(v1)*np.linalg.norm(v2))
                cos_angles.append(cos)

        q = 1 - (3/8)*np.sum((np.array(cos_angles) + 1/3)**2)
        q_frame.append(q)

    q_list.append(np.mean(q_frame) if q_frame else np.nan)

    # ---------- n_a (accepted H-bonds) ----------
    frame_hbonds = hbonds[hbonds[:,0] == ts.frame]

    acceptors = frame_hbonds[:,3].astype(int)
    counts = np.bincount(acceptors, minlength=len(O))

    na_list.append(np.mean(counts))

    # ---------- build adjacency graph ----------
    graph = {i: [] for i in range(len(O))}

    for hb in frame_hbonds:
        donor = int(hb[1])
        acceptor = int(hb[3])
        graph[donor].append(acceptor)

    # ---------- w4 (chain length of 4 waters) ----------
    lengths = []

    def dfs(path):
        if len(path) == 4:
            # compute total length
            coords = O_pos[path]
            dist = np.sum(np.linalg.norm(np.diff(coords, axis=0), axis=1))
            lengths.append(dist)
            return

        last = path[-1]
        for nxt in graph[last]:
            if nxt not in path:
                dfs(path + [nxt])

    for start in graph:
        dfs([start])

    w4_list.append(np.mean(lengths) if lengths else np.nan)

    # ---------- qcos (alignment along wire) ----------
    angles = []

    for path_len4 in lengths[:50]:  # limit for speed
        pass

    # simpler proxy: angle between nearest neighbors
    for i, neigh in enumerate(neighbors):
        if len(neigh) < 3:
            continue

        v1 = O_pos[neigh[1]] - O_pos[i]
        v2 = O_pos[neigh[2]] - O_pos[i]

        cos = np.dot(v1, v2) / (np.linalg.norm(v1)*np.linalg.norm(v2))
        angles.append(cos)

    qcos_list.append(np.mean(angles) if angles else np.nan)

# =========================
# CONVERT TO ARRAYS
# =========================
q_arr = np.array(q_list)
na_arr = np.array(na_list)
w4_arr = np.array(w4_list)
qcos_arr = np.array(qcos_list)

time = np.arange(len(q_arr)) * u.trajectory.dt

# =========================
# SAVE DATA
# =========================
data = np.column_stack((time, q_arr, na_arr, w4_arr, qcos_arr))

np.savetxt(
    f"{output_prefix}_metrics.dat",
    data,
    header="time(ps)  q  n_a  w4(Å)  qcos",
    fmt="%.6f"
)

# =========================
# SUMMARY FILE
# =========================
with open(f"{output_prefix}_summary.txt", "w") as f:
    f.write("==== SUMMARY ====\n")
    f.write(f"q mean: {np.nanmean(q_arr):.3f}\n")
    f.write(f"n_a mean: {np.nanmean(na_arr):.3f}\n")
    f.write(f"w4 mean: {np.nanmean(w4_arr):.3f} Å\n")
    f.write(f"qcos mean: {np.nanmean(qcos_arr):.3f}\n")

# =========================
# PLOTS
# =========================
plt.figure()
plt.plot(time, q_arr)
plt.xlabel("Time (ps)")
plt.ylabel("q (tetrahedral)")
plt.savefig(f"{output_prefix}_q.png")

plt.figure()
plt.plot(time, na_arr)
plt.xlabel("Time (ps)")
plt.ylabel("n_a")
plt.savefig(f"{output_prefix}_na.png")

plt.figure()