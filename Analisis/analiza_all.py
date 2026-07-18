"""
TIP3P water analysis at 277 K (NVT)
Computes:
  1. Tetrahedral order parameter q_tet
  2. Self-diffusion coefficient D (Einstein MSD)
  3. Hydrogen bond lifetime — intermittent C(t)

Input:  water_tip3p.pdb  +  tip3p_277K_NVT.dcd
Output: tip3p_results.png  +  printed summary
"""

import warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from scipy import stats
from scipy.ndimage import gaussian_filter1d
import MDAnalysis as mda
from MDAnalysis.analysis.msd import EinsteinMSD
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.lib.distances import distance_array, calc_angles, calc_bonds
warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────
PDB = "water_tip5p.pdb"
DCD = "tip5p_295K_NPT.dcd"
OUT = "anal/tip5p_295K_results.png"

# ── H-bond criteria (Luzar-Chandler) ──────────────────────────────────
HBOND_OO   = 3.5    # Å  donor O ... acceptor O
HBOND_ANG  = 150.0  # deg  D-H...A minimum angle
STRIDE_HB  = 5      # analyse every 5th frame for HBA
STRIDE_Q   = 5      # every 5th frame for q_tet
MAX_TAU    = 200    # max lag frames for ACF

# MSD fit window
MSD_LO, MSD_HI = 0.10, 0.90

# ── Colours ────────────────────────────────────────────────────────────
COLOR  = "#4FC3F7"
bg, axc, tc = "#0D1117", "#161B22", "#C9D1D9"

def style_ax(ax, title, xlabel, ylabel):
    ax.set_facecolor(axc)
    ax.tick_params(colors=tc, labelsize=9)
    ax.xaxis.label.set_color(tc); ax.yaxis.label.set_color(tc)
    ax.title.set_color("#E6EDF3")
    for sp in ax.spines.values(): sp.set_color("#30363D")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)

# ══════════════════════════════════════════════════════════════════════
print("Loading trajectory...")
u = mda.Universe(PDB, DCD)
n_frames = u.trajectory.n_frames
dt_ps    = u.trajectory.dt          # 10 ps
n_water  = u.atoms.n_residues
box0     = u.trajectory[0].dimensions
print(f"  {n_water} water molecules | {n_frames} frames | dt={dt_ps:.1f} ps")
print(f"  Box: {box0[:3]} Å")
print(f"  Total time: {n_frames * dt_ps / 1000:.1f} ns")

# ══════════════════════════════════════════════════════════════════════
# 1. TETRAHEDRAL ORDER PARAMETER
# ══════════════════════════════════════════════════════════════════════
print("\n[1/3] Tetrahedral order parameter q_tet ...")

def compute_qtet(u, stride=STRIDE_Q, n_neigh=4, cutoff=4.0):
    ow    = u.select_atoms("name O")
    q_arr = []
    t_arr = []
    frames = u.trajectory[::stride]
    for k, ts in enumerate(frames):
        pos  = ow.positions.copy()
        box  = ts.dimensions
        box3 = box[:3]
        dmat = distance_array(pos, pos, box=box)
        q_frame = []
        for i in range(len(pos)):
            row = dmat[i].copy(); row[i] = np.inf
            nn_idx = np.argpartition(row, n_neigh)[:n_neigh]
            if row[nn_idx].max() > cutoff:
                continue
            vecs = []
            for j in nn_idx:
                dr   = pos[j] - pos[i]
                dr  -= np.round(dr / box3) * box3
                norm = np.linalg.norm(dr)
                if norm > 0:
                    vecs.append(dr / norm)
            if len(vecs) < n_neigh:
                continue
            q_i = 0.0
            for a in range(n_neigh):
                for b in range(a+1, n_neigh):
                    cos_ab = np.clip(np.dot(vecs[a], vecs[b]), -1, 1)
                    q_i   += (cos_ab + 1.0/3.0)**2
            q_frame.append(1.0 - (3.0/8.0) * q_i)
        q_arr.append(np.mean(q_frame))
        t_arr.append(ts.time)
        if k % 20 == 0:
            print(f"  frame {k}/{len(frames)}  q={q_arr[-1]:.4f}")
    u.trajectory.rewind()
    return np.array(t_arr), np.array(q_arr)

t_q, q_arr = compute_qtet(u)
q_mean = q_arr.mean()
q_std  = q_arr.std()
print(f"  <q_tet> = {q_mean:.4f} ± {q_std:.4f}")

# ══════════════════════════════════════════════════════════════════════
# 2. SELF-DIFFUSION COEFFICIENT (MSD)
# ══════════════════════════════════════════════════════════════════════
print("\n[2/3] Self-diffusion coefficient (Einstein MSD) ...")

# Need nojump trajectory — build it from scratch using unwrapped positions
# (DCD may have wrapped coords; unwrap manually)
ow = u.select_atoms("name O")
n_ow = ow.n_atoms

positions_all = np.zeros((n_frames, n_ow, 3))
box_all       = np.zeros((n_frames, 3))

print("  Reading positions...")
for k, ts in enumerate(u.trajectory):
    positions_all[k] = ow.positions.copy()
    box_all[k]       = ts.dimensions[:3]
u.trajectory.rewind()

# Unwrap: remove PBC jumps frame-by-frame
print("  Unwrapping PBC...")
pos_unwrap = positions_all.copy()
for k in range(1, n_frames):
    dr    = pos_unwrap[k] - pos_unwrap[k-1]
    box_k = box_all[k]
    dr   -= np.round(dr / box_k) * box_k
    pos_unwrap[k] = pos_unwrap[k-1] + dr

# Compute MSD
print("  Computing MSD...")
t_msd = np.arange(n_frames) * dt_ps
msd   = np.zeros(n_frames)
for tau in range(n_frames):
    if tau % 200 == 0:
        print(f"  tau {tau}/{n_frames}")
    disp  = pos_unwrap[tau:] - pos_unwrap[:n_frames-tau]
    msd[tau] = np.mean(np.sum(disp**2, axis=2))

# Linear fit on middle portion
lo  = int(n_frames * MSD_LO)
hi  = int(n_frames * MSD_HI)
slope, intercept, r_val, _, _ = stats.linregress(t_msd[lo:hi], msd[lo:hi])
D_m2s = (slope / 6.0) * 1e-8   # Å²/ps → m²/s
print(f"  D = {D_m2s:.3e} m²/s  (R²={r_val**2:.5f})")

# ══════════════════════════════════════════════════════════════════════
# 3. HYDROGEN BOND LIFETIME — INTERMITTENT C(t)
# ══════════════════════════════════════════════════════════════════════
print("\n[3/3] Hydrogen bond lifetime (intermittent ACF) ...")

hba = HydrogenBondAnalysis(
    universe=u,
    donors_sel="name O",
    hydrogens_sel="name H1 H2",
    acceptors_sel="name O",
    d_h_cutoff=1.2,
    d_a_cutoff=HBOND_OO,
    d_h_a_angle_cutoff=HBOND_ANG,
    update_selections=False,
)
hba.run(step=STRIDE_HB, verbose=True)
u.trajectory.rewind()

hb      = hba.results.hbonds
print(f"  Total H-bond instances: {len(hb)}")

# Manual intermittent ACF
frames_hb  = hb[:, 0].astype(int)
donors_hb  = hb[:, 1].astype(int)
acceptors_hb = hb[:, 3].astype(int)

uniq_f  = np.unique(frames_hb)
n_f     = len(uniq_f)
f2i     = {f: i for i, f in enumerate(uniq_f)}

pairs = {}
for fr, d, a in zip(frames_hb, donors_hb, acceptors_hb):
    key = (int(d), int(a))
    pairs.setdefault(key, set()).add(f2i[int(fr)])

tau_max = min(MAX_TAU, n_f // 2)
num = np.zeros(tau_max)
den = 0.0
for key, ps in pairs.items():
    h = np.zeros(n_f, dtype=np.float32)
    for idx in ps:
        h[idx] = 1.0
    if h[0] == 0:
        continue
    den += 1.0
    for tau in range(tau_max):
        num[tau] += h[tau]

acf    = num / den
acf   /= acf[0]
dt_hb  = dt_ps * STRIDE_HB
tau_ps = np.arange(tau_max) * dt_hb
hb_life = float(trapezoid(acf, tau_ps))
print(f"  τ_HB (intermittent) = {hb_life:.2f} ps")

# ══════════════════════════════════════════════════════════════════════
# PLOTTING — 3-panel figure
# ══════════════════════════════════════════════════════════════════════
print("\nPlotting...")

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
fig.patch.set_facecolor(bg)

# ── Panel 1: q_tet over time ──────────────────────────────────────────
ax = axes[0]
style_ax(ax, "Tetrahedral Order Parameter", "t (ns)", "q_tet")
t_ns = t_q / 1000
win  = max(1, len(q_arr) // 20)
ax.plot(t_ns, q_arr, color=COLOR, lw=0.8, alpha=0.35)
q_smooth = np.convolve(q_arr, np.ones(win)/win, mode='valid')
t_smooth = t_ns[win-1:]
ax.plot(t_smooth, q_smooth, color=COLOR, lw=2.2)
ax.axhline(q_mean, color="#FFC107", lw=1.4, ls="--",
           label=f"⟨q⟩ = {q_mean:.4f} ± {q_std:.4f}")
ax.set_ylim(0, 1)
ax.legend(framealpha=0.2, facecolor="#21262D", labelcolor=tc, fontsize=9)

# ── Panel 2: MSD + linear fit ─────────────────────────────────────────
ax = axes[1]
style_ax(ax, "Mean Squared Displacement", "t (ps)", "MSD (Å²)")
ax.plot(t_msd, msd, color=COLOR, lw=1.5, alpha=0.7, label="MSD")
t_fit = t_msd[lo:hi]
ax.plot(t_fit, slope * t_fit + intercept,
        color="#FF7043", lw=2.5, ls="--",
        label=f"Linear fit (R²={r_val**2:.4f})\nD = {D_m2s:.3e} m²/s")
ax.legend(framealpha=0.2, facecolor="#21262D", labelcolor=tc, fontsize=9)

# ── Panel 3: H-bond ACF ───────────────────────────────────────────────
ax = axes[2]
style_ax(ax, "Intermittent H-bond ACF C(t)", "t (ps)", "C(t)")
ax.plot(tau_ps, acf, color=COLOR, lw=2.0,
        label=f"TIP3P 277 K\nτ = {hb_life:.1f} ps")
ax.set_ylim(0, 1); ax.set_xlim(0, tau_ps[-1])
ax.axhline(1/np.e, color="#484F58", lw=0.8, ls=":",
           label="1/e ≈ 0.368")
ax.legend(framealpha=0.2, facecolor="#21262D", labelcolor=tc, fontsize=9)

fig.suptitle("TIP3P Water — 277 K NVT Analysis",
             color="#E6EDF3", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(OUT, dpi=180, bbox_inches="tight", facecolor=bg)
print(f"Saved → {OUT}")

# ── Summary ────────────────────────────────────────────────────────────
print("\n" + "="*50)
print("  SUMMARY ")
print("="*50)
print(f"  q_tet          = {q_mean:.4f} ± {q_std:.4f}")
print(f"  D              = {D_m2s:.3e} m²/s")
print(f"  τ_HB (interm.) = {hb_life:.2f} ps")
print(f"  n_water        = {n_water}")
print(f"  T (log mean)   = ~277 K (NVT target)")
print("="*50)