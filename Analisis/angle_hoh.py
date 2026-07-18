import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

# =========================
# INPUT
# =========================
u = mda.Universe("water_tip3p.pdb", "amoeba_260K_NPT.dcd")
output_file = "HOH_angles_amoeba_260K_NPT.dat"
# =========================
# SELECT WATER
# =========================
waters = u.select_atoms("resname HOH or resname WAT")

# =========================
# COLLECT ANGLES
# =========================
angles = []

for ts in u.trajectory:
    for res in waters.residues:
        O = res.atoms.select_atoms("name O")
        H = res.atoms.select_atoms("name H1 H2 H")

        if len(O) != 1 or len(H) != 2:
            continue  # safety

        O = O.positions[0]
        H1, H2 = H.positions

        v1 = H1 - O
        v2 = H2 - O

        cos_theta = np.dot(v1, v2) / (
            np.linalg.norm(v1) * np.linalg.norm(v2)
        )

        # numerical safety
        cos_theta = np.clip(cos_theta, -1.0, 1.0)

        theta = np.degrees(np.arccos(cos_theta))
        angles.append(theta)

angles = np.array(angles)

# =========================
# SAVE DATA
# =========================
np.savetxt(
    output_file,
    angles,
    header="HOH angle [deg]",
    comments=""
)

print(f"Zapisano {len(angles)} kątów HOH")

print("Min angle:", np.min(angles))
print("Max angle:", np.max(angles))
print("Mean angle:", np.mean(angles))
print("Std angle:", np.std(angles))
print("NaNs:", np.isnan(angles).sum())
