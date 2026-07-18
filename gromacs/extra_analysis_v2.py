import os, warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid as trapz
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.lib.distances import calc_angles, calc_bonds
warnings.filterwarnings("ignore")

MODELS = {
    "FBA/eps":     {"tpr": "fbaeps/npt_correct/npt.tpr",  "xtc": "fbaeps/npt_correct/npt.xtc"},
    "SPC/Fw":      {"tpr": "spcfw/npt/npt_fixed.tpr",     "xtc": "spcfw/npt/npt.xtc"},
    "TIP4P/2005f": {"tpr": "tip4p2005f/npt/npt.tpr",      "xtc": "tip4p2005f/npt/npt.xtc"},
}
COLORS  = {"FBA/eps": "#2196F3", "SPC/Fw": "#E91E63", "TIP4P/2005f": "#4CAF50"}
OUT_DIR = "water_analysis_results"
os.makedirs(OUT_DIR, exist_ok=True)

HBOND_OO_DIST   = 3.5    # Å  O...O cutoff
HBOND_ANGLE     = 150.0  # deg D-H...A minimum
STRIDE_HBOND    = 20
STRIDE_GEOM     = 10
MAX_TAU_FRAMES  = 150    # max lag for ACF

bg, axc, tc = "#0D1117", "#161B22", "#C9D1D9"

def style_ax(ax):
    ax.set_facecolor(axc)
    ax.tick_params(colors=tc, labelsize=8)
    ax.xaxis.label.set_color(tc)
    ax.yaxis.label.set_color(tc)
    ax.title.set_color("#E6EDF3")
    for sp in ax.spines.values():
        sp.set_color("#30363D")

# ── FIX 1: RDF with correct ranges to skip intramolecular peaks ───────
def read_rdf(u, sel1, sel2, nbins=300, rng=(1.5, 6.0), excl=None):
    ag1 = u.select_atoms(sel1)
    ag2 = u.select_atoms(sel2)
    kwargs = dict(nbins=nbins, range=rng)
    if excl:
        kwargs["exclusion_block"] = excl
    rdf = InterRDF(ag1, ag2, **kwargs)
    rdf.run(verbose=False)
    u.trajectory.rewind()
    return rdf.results.bins, rdf.results.rdf

def get_oh_hoh_geometry(u):
    oh_lengths, hoh_angles = [], []
    traj = u.trajectory[::STRIDE_GEOM]
    print(f"    {len(traj)} frames")
    for ts in traj:
        box = ts.dimensions
        for mol in u.residues:
            ow = mol.atoms.select_atoms("name OW")
            hw = mol.atoms.select_atoms("name HW*")
            if len(ow) != 1 or len(hw) < 2:
                continue
            o, h1, h2 = ow.positions[0], hw.positions[0], hw.positions[1]
            oh_lengths.extend([
                calc_bonds(o[None], h1[None], box=box)[0],
                calc_bonds(o[None], h2[None], box=box)[0],
            ])
            hoh_angles.append(np.degrees(
                calc_angles(h1[None], o[None], h2[None], box=box)[0]))
    u.trajectory.rewind()
    return np.array(oh_lengths), np.array(hoh_angles)

def run_hba(u):
    hba = HydrogenBondAnalysis(
        universe=u,
        donors_sel="name OW",
        hydrogens_sel="name HW*",
        acceptors_sel="name OW",
        d_h_cutoff=1.2,
        d_a_cutoff=HBOND_OO_DIST,
        d_h_a_angle_cutoff=HBOND_ANGLE,
        update_selections=False,
    )
    hba.run(step=STRIDE_HBOND, verbose=True)
    u.trajectory.rewind()
    return hba

# ── FIX 2: Debug and correctly extract H-bond distances ───────────────
def extract_hbond_geometry(hba, u):
    """
    HBA results columns: frame, donor_idx, H_idx, acceptor_idx, dist, angle
    dist = H...acceptor distance (Å),  angle = D-H...A (degrees)
    """
    hbonds = hba.results.hbonds
    if len(hbonds) == 0:
        return np.array([]), np.array([])

    h_a_dist = hbonds[:, 4].astype(float)   # H...O acceptor distance
    dha_angle = hbonds[:, 5].astype(float)  # D-H...A angle

    # Sanity print
    print(f"    H...O dist: min={h_a_dist.min():.3f}  max={h_a_dist.max():.3f}  "
          f"mean={h_a_dist.mean():.3f} A")
    print(f"    D-H...A angle: min={dha_angle.min():.1f}  max={dha_angle.max():.1f} deg")
    return h_a_dist, dha_angle

# ── FIX 3: Manual intermittent ACF (no hba.lifetime() dependency) ─────
def compute_intermittent_acf(hba, u, max_tau=MAX_TAU_FRAMES):
    """
    Intermittent C(t): probability that a pair bonded at t=0
    is STILL bonded at t (allowing transient breaks).
    Vectorized over all observed pairs.
    """
    hbonds  = hba.results.hbonds
    if len(hbonds) == 0:
        return np.array([]), np.array([])

    frames_all = hbonds[:, 0].astype(int)
    donors_all = hbonds[:, 1].astype(int)
    acceptors_all = hbonds[:, 3].astype(int)

    uniq_frames = np.unique(frames_all)
    n_frames = len(uniq_frames)
    frame_to_idx = {f: i for i, f in enumerate(uniq_frames)}

    # For each unique pair build presence array h[frame_idx]
    pairs = {}
    for fnum, d, a in zip(frames_all, donors_all, acceptors_all):
        key = (int(d), int(a))
        fidx = frame_to_idx[int(fnum)]
        if key not in pairs:
            pairs[key] = set()
        pairs[key].add(fidx)

    max_tau = min(max_tau, n_frames // 2)
    acf_num = np.zeros(max_tau)
    acf_den = 0.0

    for key, present_set in pairs.items():
        h = np.zeros(n_frames, dtype=np.float32)
        for idx in present_set:
            h[idx] = 1.0
        if h[0] == 0:
            continue
        acf_den += 1.0
        # Intermittent: h(0)*h(t) for each origin t0=0 per pair
        for tau in range(max_tau):
            acf_num[tau] += h[tau]

    if acf_den == 0:
        return np.array([]), np.array([])

    acf = acf_num / acf_den
    acf = acf / acf[0]   # normalise so C(0)=1

    # Time axis
    dt_ps = u.trajectory.dt * STRIDE_HBOND
    tau_ps = np.arange(max_tau) * dt_ps

    hb_life = float(trapz(acf, tau_ps))
    print(f"    tau_HB (intermittent) = {hb_life:.2f} ps")
    return tau_ps, acf

# ══════════════════════════════════════════════════════════════════════
results = {}
for label, paths in MODELS.items():
    print(f"\n{'='*50}\n  {label}\n{'='*50}")
    u = mda.Universe(paths["tpr"], paths["xtc"])

    print("[1/4] RDF O-H (start 1.5 A) and H-H (start 2.0 A) ...")
    # Skip intramolecular peak by starting range beyond bond length
    r_oh, g_oh = read_rdf(u, "name OW", "name HW*", rng=(1.5, 5.0))
    r_hh, g_hh = read_rdf(u, "name HW*", "name HW*", rng=(2.0, 5.0))

    print("[2/4] Intramolecular geometry ...")
    oh_len, hoh_ang = get_oh_hoh_geometry(u)
    print(f"    OH mean={oh_len.mean():.4f} A  HOH mean={hoh_ang.mean():.2f} deg")

    print("[3/4] H-bond analysis ...")
    hba = run_hba(u)
    h_oa_dist, dha_ang = extract_hbond_geometry(hba, u)

    print("[4/4] Intermittent ACF (manual) ...")
    tau_ps, acf = compute_intermittent_acf(hba, u)
    hb_life = float(trapz(acf, tau_ps)) if len(acf) > 1 else float('nan')

    results[label] = {
        "rdf_oh":    (r_oh, g_oh),
        "rdf_hh":    (r_hh, g_hh),
        "oh_len":    oh_len,
        "hoh_ang":   hoh_ang,
        "h_oa_dist": h_oa_dist,
        "dha_angle": dha_ang,
        "hb_acf":    (tau_ps, acf),
        "hb_life":   hb_life,
    }

# ── Figure 1: RDF O-H and H-H ─────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.patch.set_facecolor(bg)
for ax, key, title, xl, rng in zip(
        axes,
        ["rdf_oh", "rdf_hh"],
        ["O-H Radial Distribution Function (intermolecular)",
         "H-H Radial Distribution Function (intermolecular)"],
        ["r (A)", "r (A)"],
        [(1.5, 5.0), (2.0, 5.0)]):
    style_ax(ax)
    for label, res in results.items():
        r, g = res[key]
        ax.plot(r, g, color=COLORS[label], lw=1.8, label=label)
    ax.axhline(1.0, color="#484F58", lw=0.8, ls="--")
    ax.set_xlabel(xl); ax.set_ylabel("g(r)")
    ax.set_title(title)
    ax.set_xlim(rng)
    ax.legend(framealpha=0.2, facecolor="#21262D", labelcolor=tc, fontsize=8)
fig.suptitle("Partial RDFs — Intermolecular", color="#E6EDF3",
             fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/rdf_OH_HH.png", dpi=180, bbox_inches="tight", facecolor=bg)
print("\nSaved rdf_OH_HH.png")

# ── Figure 2: Intramolecular geometry ─────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.patch.set_facecolor(bg)
for ax in axes: style_ax(ax)
for label, res in results.items():
    c = COLORS[label]
    d = res["oh_len"]
    h, e = np.histogram(d, bins=np.linspace(0.85, 1.15, 150), density=True)
    axes[0].plot(0.5*(e[:-1]+e[1:]), h, color=c, lw=1.8,
                 label=f"{label}  mean={d.mean():.4f} A  std={d.std():.4f}")
    a = res["hoh_ang"]
    h, e = np.histogram(a, bins=np.linspace(90, 130, 120), density=True)
    axes[1].plot(0.5*(e[:-1]+e[1:]), h, color=c, lw=1.8,
                 label=f"{label}  mean={a.mean():.2f}  std={a.std():.2f} deg")
axes[0].set_xlabel("O-H bond length (A)"); axes[0].set_ylabel("Probability density")
axes[0].set_title("O-H Bond Length Distribution")
axes[0].legend(framealpha=0.2, facecolor="#21262D", labelcolor=tc, fontsize=8)
axes[1].set_xlabel("H-O-H angle (deg)"); axes[1].set_ylabel("Probability density")
axes[1].set_title("H-O-H Angle Distribution")
axes[1].legend(framealpha=0.2, facecolor="#21262D", labelcolor=tc, fontsize=8)
fig.suptitle("Intramolecular Geometry", color="#E6EDF3", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/intramolecular_geometry.png", dpi=180,
            bbox_inches="tight", facecolor=bg)
print("Saved intramolecular_geometry.png")

# ── Figure 3: H-bond geometry ─────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.patch.set_facecolor(bg)
for ax in axes: style_ax(ax)
for label, res in results.items():
    c = COLORS[label]
    d = res["h_oa_dist"]
    a = res["dha_angle"]
    if len(d) == 0:
        print(f"  No H-bonds found for {label}")
        continue
    h, e = np.histogram(d, bins=np.linspace(1.4, 2.8, 120), density=True)
    axes[0].plot(0.5*(e[:-1]+e[1:]), h, color=c, lw=1.8,
                 label=f"{label}  mean={d.mean():.3f} A")
    h, e = np.histogram(a, bins=np.linspace(140, 180, 100), density=True)
    axes[1].plot(0.5*(e[:-1]+e[1:]), h, color=c, lw=1.8,
                 label=f"{label}  mean={a.mean():.1f} deg")
axes[0].set_xlabel("H...O distance (A)"); axes[0].set_ylabel("Probability density")
axes[0].set_title("H-bond H...O Distance Distribution")
axes[0].legend(framealpha=0.2, facecolor="#21262D", labelcolor=tc, fontsize=8)
axes[1].set_xlabel("D-H...A angle (deg)"); axes[1].set_ylabel("Probability density")
axes[1].set_title("H-bond D-H...A Angle Distribution")
axes[1].legend(framealpha=0.2, facecolor="#21262D", labelcolor=tc, fontsize=8)
fig.suptitle("Hydrogen Bond Geometry", color="#E6EDF3", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/hbond_geometry.png", dpi=180, bbox_inches="tight", facecolor=bg)
print("Saved hbond_geometry.png")

# ── Figure 4: H-bond ACF ──────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 5))
fig.patch.set_facecolor(bg); style_ax(ax)
any_plotted = False
for label, res in results.items():
    tau_ps, acf = res["hb_acf"]
    lt = res["hb_life"]
    if len(tau_ps) < 2:
        print(f"  No ACF data for {label}")
        continue
    ax.plot(tau_ps, acf, color=COLORS[label], lw=2.0,
            label=f"{label}  tau={lt:.1f} ps")
    any_plotted = True
if not any_plotted:
    ax.text(0.5, 0.5, "ACF data empty — check HBA output",
            transform=ax.transAxes, color=tc, ha="center")
ax.set_xlabel("t (ps)"); ax.set_ylabel("C(t)")
ax.set_title("Intermittent H-bond Autocorrelation C(t)")
ax.set_ylim(0, 1)
ax.legend(framealpha=0.2, facecolor="#21262D", labelcolor=tc, fontsize=9)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/hbond_acf.png", dpi=180, bbox_inches="tight", facecolor=bg)
print("Saved hbond_acf.png")

print("\nAll done.")
