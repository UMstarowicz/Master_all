import numpy as np
import MDAnalysis as mda
import os, shutil

GRO      = 'gromacs_input/conf.gro'
TRR      = 'gromacs_input/sample_flex.trr'
LOAD_DIR = 'load'

INTERFACES = [0.103, 0.105, 0.107, 0.108, 0.109,
              0.110, 0.111, 0.112, 0.113, 0.115]
TARGETS = [INTERFACES[0], INTERFACES[0]] + list(INTERFACES[1:])

def find_lambda(ow_pos, hw_pos, box):
    max_oh = 0.0
    for i in range(len(ow_pos)):
        for h_local in range(2):
            h = hw_pos[2*i + h_local].copy()
            d = h - ow_pos[i]
            for dim in range(3):
                L = box[dim]
                d[dim] -= L * round(d[dim] / L)
            r = float(np.linalg.norm(d))
            if r > max_oh:
                max_oh = r
    return max_oh

print("Loading trajectory...")
u  = mda.Universe(GRO, TRR)
OW = u.select_atoms('name OW')
HW = u.select_atoms('name HW1 or name HW2')

lambdas = []
for ts in u.trajectory:
    ow  = OW.positions / 10.0
    hw  = HW.positions / 10.0
    box = ts.dimensions[:3] / 10.0
    lv  = find_lambda(ow, hw, box)
    lambdas.append((ts.frame, lv))

lambdas = np.array(lambdas)
print(f"Scanned {len(lambdas)} frames")
print(f"λ range: {lambdas[:,1].min():.4f} – {lambdas[:,1].max():.4f} nm")

for ens_idx, target in enumerate(TARGETS):
    folder      = os.path.join(LOAD_DIR, f'{ens_idx:03d}')
    accepted    = os.path.join(folder, 'accepted')
    os.makedirs(accepted, exist_ok=True)

    diffs       = np.abs(lambdas[:, 1] - target)
    best_frame  = int(lambdas[np.argmin(diffs), 0])
    best_lambda = float(lambdas[np.argmin(diffs), 1])

    # 1. Write conf.gro
    u.trajectory[best_frame]
    gro_path = os.path.join(folder, 'conf.gro')
    with mda.Writer(gro_path, u.select_atoms('all').n_atoms) as W:
        W.write(u.select_atoms('all'))

    # 2. Copy conf.gro to accepted/0.000
    traj_frame = os.path.join(accepted, '0.000')
    shutil.copy(gro_path, traj_frame)

    # 3. Write order.txt
    with open(os.path.join(folder, 'order.txt'), 'w') as f:
        f.write('# Re-calculated order parameters.\n')
        f.write('#     Time       Orderp\n')
        f.write(f'         0     {best_lambda:.6f}\n')

    # 4. Write traj.txt
    with open(os.path.join(folder, 'traj.txt'), 'w') as f:
        f.write('# Trajectory file\n')
        f.write(f'0  0.000  accepted/0.000\n')

    print(f"  Ensemble {ens_idx:03d}: target={target:.3f} nm  "
          f"→ frame {best_frame}, λ={best_lambda:.4f} nm")

print("\nDone.")
