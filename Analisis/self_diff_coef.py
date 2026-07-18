import MDAnalysis as mda
import numpy as np

u = mda.Universe('top.psf','traj.dcd')
oxygen = u.select_atoms('resname TIP3 and name O')
msd = np.zeros(len(u.trajectory))
count = 0
for ts in u.trajectory:  # tutaj u.trajectory.all_frames() może być duże; lepiej wybrać co kilka klatek
    # opcjonalnie unwrap: ts.unwrap()
    dr = oxygen.positions - oxygen.positions  # na początku t0
    # wyznaczanie dr dla kolejnych tau
    # W MDAnalysis 2.0: można użyć DiffusionAnalysis
