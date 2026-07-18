import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import distances

# =====================================
# LOAD TRAJECTORY
# =====================================

u = mda.Universe("prod.tpr", "prod.trr")

oxygens = u.select_atoms("name OW")

q_values = []

# =====================================
# LOOP OVER FRAMES
# =====================================

for ts in u.trajectory[5000::10]:

    positions = oxygens.positions

    dist_matrix = distances.distance_array(
        positions,
        positions,
        box=u.dimensions
    )

    np.fill_diagonal(dist_matrix, np.inf)

    for i in range(len(oxygens)):

        # four nearest neighbours
        nn_idx = np.argsort(dist_matrix[i])[:4]

        neigh = positions[nn_idx]
        center = positions[i]

        vectors = neigh - center

        # normalize
        vectors /= np.linalg.norm(vectors, axis=1)[:, None]

        q = 0.0

        for j in range(3):
            for k in range(j+1, 4):

                cospsi = np.dot(vectors[j], vectors[k])

                q += (cospsi + 1/3)**2

        q = 1 - (3/8)*q

        q_values.append(q)

# =====================================
# RESULTS
# =====================================
q_values = np.array(q_values)
q_mean = np.mean(q_values)

print(f"Average tetrahedral order parameter q = {q_mean:.4f}")
print("q min/max check:", np.min(q_values), np.max(q_values))
