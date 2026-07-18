import MDAnalysis as mda
import numpy as np
from scipy.stats import linregress

u = mda.Universe('water_tip5p.pdb', 'tip5p_260K_NPT.dcd')

water_O = u.select_atoms('element O')
print(len(water_O))
tetra_values = []
positions = []
hb_series = []

dt = u.trajectory.dt  # ps (OpenMM provides this sometimes)

def tetrahedral_order(positions, box):
    N = len(positions)
    q_vals = []

    for i in range(N):
        rij = positions - positions[i]
        rij -= box * np.round(rij / box)

        dist = np.linalg.norm(rij, axis=1)

        dist[i] = np.inf

        nn = np.argsort(dist)[:4]

        if len(nn) < 4:
            q_vals.append(np.nan)
            continue

        vecs = rij[nn]
        vecs /= np.linalg.norm(vecs, axis=1)[:, None]

        q = 1.0

        for a in range(3):
            for b in range(a + 1, 4):
                cos = np.dot(vecs[a], vecs[b])
                q -= (3/8) * (cos + 1/3)**2

        q_vals.append(q)

    return np.array(q_vals)

def tetrahedral_order(positions, box):
    N = len(positions)
    q_vals = []

    for i in range(N):
        rij = positions - positions[i]

        # PBC correction (MDAnalysis-friendly)
        rij -= box * np.round(rij / box)

        dist = np.linalg.norm(rij, axis=1)
        dist[i] = np.inf

        nn = np.argsort(dist)[1:5]
        if len(nn) != 4:
            continue
        vecs = rij[nn]

        vecs /= np.linalg.norm(vecs, axis=1)[:, None]

        q = 1.0
        for a in range(3):
            for b in range(a + 1, 4):
                cos = np.dot(vecs[a], vecs[b])
                q -= (3/8) * (cos + 1/3)**2

        q_vals.append(q)

    return np.array(q_vals)

def hydrogen_bonds(frame_O, box, cutoff=3.5):
    N = len(frame_O)
    hb = []

    for i in range(N):
        rij = frame_O - frame_O[i]
        rij -= box * np.round(rij / box)

        dist = np.linalg.norm(rij, axis=1)

        bonded = np.where((dist < cutoff) & (dist > 0))[0]
        hb.append(set(bonded))

    return hb

for ts in u.trajectory:
    O = water_O.positions.copy()
    box = u.dimensions[:3]

    positions.append(O.copy())

    q_frame = tetrahedral_order(O, box)
    if len(q_frame) > 0:
        tetra_values.append(np.nanmean(q_frame))
    else:
        tetra_values.append(np.nan)

    hb_series.append(hydrogen_bonds(O, box))

def hydrogen_bond_matrix(frame_O, box, cutoff=3.5):
    rij = frame_O[:, None, :] - frame_O[None, :, :]
    rij -= box * np.round(rij / box)

    dist = np.linalg.norm(rij, axis=-1)

    hb = (dist < cutoff) & (dist > 0)
    return hb.astype(int)

def hb_lifetime(hb_series, max_t):
    N = len(hb_series)
    C = np.zeros(max_t)

    for t in range(max_t):
        num = 0
        den = 0

        for i in range(N - t):
            A = hb_series[i]
            B = hb_series[i + t]

            num += np.sum(A * B)
            den += np.sum(A)

        C[t] = num / den if den > 0 else 0

    return C

u.trajectory.add_transformations(
    mda.transformations.unwrap(u.atoms)
)

def coordination_number(frame_O, box, cutoff=3.5):
    rij = frame_O[:, None, :] - frame_O[None, :, :]
    rij -= box * np.round(rij / box)

    dist = np.linalg.norm(rij, axis=-1)

    np.fill_diagonal(dist, np.inf)

    cn = np.sum(dist < cutoff, axis=1)

    return cn

positions = []

for ts in u.trajectory:
    positions.append(water_O.positions.copy())

def msd(positions):
    n_frames = len(positions)

    r0 = positions[0].copy()
    msd_vals = np.zeros(n_frames)

    for t in range(n_frames):
        dr = positions[t] - r0
        msd_vals[t] = np.mean(np.sum(dr**2, axis=1))

    return msd_vals

box = u.dimensions[:3]

msd_vals = msd(positions)

t = np.arange(len(msd_vals)) * dt

large_t = int(len(msd_vals) * 0.3)

mask = np.isfinite(msd_vals)
t_clean = t[mask]
msd_clean = msd_vals[mask]

slope, _, _, _, _ = linregress(t_clean[large_t:], msd_clean[large_t:])

D = slope / 6
print("Diffusion coefficient:", D)

