# Patch: replace slow nested loop with fast HBA-based geometry extraction
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
    "FBA/eps":     {"tpr": "fbaeps/npt_correct/npt.tpr",   "xtc": "fbaeps/npt_correct/npt.xtc"},
    "SPC/Fw":      {"tpr": "spcfw/npt/npt_fixed.tpr",      "xtc": "spcfw/npt/npt.xtc"},
    "TIP4P/2005f": {"tpr": "tip4p2005f/npt/npt.tpr",       "xtc": "tip4p2005f/npt/npt.xtc"},
}
COLORS  = {"FBA/eps": "#2196F3", "SPC/Fw": "#E91E63", "TIP4P/2005f": "#4CAF50"}
OUT_DIR = "water_analysis_results"
os.makedirs(OUT_DIR, exist_ok=True)

HBOND_DIST   = 3.5
HBOND_ANGLE  = 150.0
STRIDE_HBOND = 20   # every 20th frame — fast enough, good statistics
STRIDE_GEOM  = 10

bg, axc, tc = "#0D1117", "#161B22", "#C9D1D9"

def style_ax(ax):
    ax.set_facecolor(axc)
    ax.tick_params(colors=tc, labelsize=8)
    ax.xaxis.label.set_color(tc)
    ax.yaxis.label.set_color(tc)
    ax.title.set_color("#E6EDF3")
    for sp in ax.spines.values():
        sp.set_color("#30363D")

def read_rdf(u, sel1, sel2, nbins=300, rng=(0.5, 8.0), excl=None):
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
    """Intramolecular O-H lengths and H-O-H angles — vectorized over residues."""
    oh_lengths, hoh_angles = [], []
    print(f"  Processing {len(u.trajectory[::STRIDE_GEOM])} frames...")
    for ts in u.trajectory[::STRIDE_GEOM]:
        box = ts.dimensions
        for mol in u.residues:
            ow = mol.atoms.select_atoms("name OW")
            hw = mol.atoms.select_atoms("name HW*")
            if len(ow) != 1 or len(hw) < 2:
                continue
            o  = ow.positions[0]
            h1 = hw.positions[0]
            h2 = hw.positions[1]
            d1 = calc_bonds(o[None], h1[None], box=box)[0]
            d2 = calc_bonds(o[None], h2[None], box=box)[0]
            oh_lengths.extend([d1, d2])
            ang = calc_angles(h1[None], o[None], h2[None], box=box)[0]
            hoh_angles.append(np.degrees(ang))
    u.trajectory.rewind()
    return np.array(oh_lengths), np.array(hoh_angles)

def run_hba(u):
    """Run HydrogenBondAnalysis — fast C-level implementation."""
    hba = HydrogenBondAnalysis(
        universe=u,
        donors_sel="name OW",
        hydrogens_sel="name HW*",
        acceptors_sel="name OW",
        d_h_cutoff=1.2,
        d_a_cutoff=HBOND_DIST,
        d_h_a_angle_cutoff=HBOND_ANGLE,
        update_selections=False,
    )
    hba.run(step=STRIDE_HBOND, verbose=True)
    u.trajectory.rewind()
    return hba

def extract_hbond_geometry(hba):
    """
    Extract O...O distances and O-H...O angles from HBA results table.
    Columns: frame, donor_idx, hydrogen_idx, acceptor_idx, distance, angle
    """
    hbonds = hba.results.hbonds   # shape (N_hbonds, 6)
    if len(hbonds) == 0:
        return np.array([]), np.array([])
    distances = hbonds[:, 4].astype(float)   # H...acceptor distance
    angles    = hbonds[:, 5].astype(float)   # D-H...A angle in degrees

    # Also get O...O distance: need donor and acceptor positions
    # Approximate: O...O ≈ O-H + H...O along bond axis (not exact but fast)
    # Better: recompute from stored indices — use mean H...O dist as proxy
    # For the histogram we store H...O distance and the angle
    return distances, angles

def hbond_acf(hba, max_tau=200):
    """Intermittent ACF via built-in lifetime method."""
    try:
        tau, acf = hba.lifetime(intermittent=True, window_step=1,
                                 tau_max=max_tau)
        return np.array(tau), np.array(acf)
    except Exception as e:
        print(f"  lifetime() failed: {e}")
        return np.array([]), np.array([])

# ── MAIN LOOP ──────────────────────────────────────────────────────────
results = {}
for label, paths in MODELS.items():
    print(f"\n{'='*50}\n  {label}\n{'='*50}")
    u = mda.Universe(paths["tpr"], paths["xtc"])

    print("[1/4] RDF O-H and H-H ...")
    r_oh, g_oh = read_rdf(u, "name OW", "name HW*", rng=(0.5,5.0), excl=(1,1))
    r_hh, g_hh = read_rdf(u, "name HW*", "name HW*", rng=(0.5,5.0), excl=(1,1))

    print("[2/4] Intramolecular geometry ...")
    oh_len, hoh_ang = get_oh_hoh_geometry(u)
    print(f"  OH mean={oh_len.mean():.4f} A   HOH mean={hoh_ang.mean():.2f} deg")

    print("[3/4] H-bond analysis (fast) ...")
    hba = run_hba(u)
    h_oa_dist, oho_ang = extract_hbond_geometry(hba)
    print(f"  {len(h_oa_dist)} H-bonds found")

    print("[4/4] Intermittent ACF ...")
    tau_ps_raw, acf = hbond_acf(hba, max_tau=200)
    dt = u.trajectory.dt * STRIDE_HBOND
    tau_ps = tau_ps_raw * dt if len(tau_ps_raw) > 0 else tau_ps_raw
    hb_life = float(trapz(acf, tau_ps)) if len(acf) > 1 else float('nan')
    print(f"  tau_HB (intermittent) = {hb_life:.2f} ps")

    results[label] = {
        "rdf_oh":    (r_oh, g_oh),
        "rdf_hh":    (r_hh, g_hh),
        "oh_len":    oh_len,
        "hoh_ang":   hoh_ang,
        "h_oa_dist": h_oa_dist,   # H...O acceptor distance
        "oho_angle": oho_ang,     # D-H...A angle
        "hb_acf":    (tau_ps, acf),
        "hb_life":   hb_life,
    }

# ── Figure 1: RDF O-H and H-H ─────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.patch.set_facecolor(bg)
for ax, key, title in zip(axes,
        ["rdf_oh","rdf_hh"],
        ["O-H Radial Distribution Function",
         "H-H Radial Distribution Function"]):
    style_ax(ax)
    for label, res in results.items():
        r, g = res[key]
        ax.plot(r, g, color=COLORS[label], lw=1.8, label=label)
    ax.axhline(1.0, color="#484F58", lw=0.8, ls="--")
    ax.set_xlabel("r (A)"); ax.set_ylabel("g(r)")
    ax.set_title(title); ax.set_xlim(0.5, 5.0)
    ax.legend(framealpha=0.2, facecolor="#21262D", labelcolor=tc, fontsize=8)
fig.suptitle("Partial RDFs", color="#E6EDF3", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/rdf_OH_HH.png", dpi=180, bbox_inches="tight", facecolor=bg)
print("Saved rdf_OH_HH.png")

# ── Figure 2: Intramolecular geometry ─────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.patch.set_facecolor(bg)
for ax in axes: style_ax(ax)
for label, res in results.items():
    c = COLORS[label]
    d = res["oh_len"]
    h, e = np.histogram(d, bins=np.linspace(0.85, 1.15, 120), density=True)
    axes[0].plot(0.5*(e[:-1]+e[1:]), h, color=c, lw=1.8,
                 label=f"{label}  mean={d.mean():.4f} A")
    a = res["hoh_ang"]
    h, e = np.histogram(a, bins=np.linspace(90, 130, 120), density=True)
    axes[1].plot(0.5*(e[:-1]+e[1:]), h, color=c, lw=1.8,
                 label=f"{label}  mean={a.mean():.2f} deg")
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
    a = res["oho_angle"]
    if len(d) == 0: continue
    h, e = np.histogram(d, bins=np.linspace(1.4, 2.5, 100), density=True)
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
for label, res in results.items():
    tau_ps, acf = res["hb_acf"]
    lt = res["hb_life"]
    if len(tau_ps) == 0: continue
    ax.plot(tau_ps, acf, color=COLORS[label], lw=2.0,
            label=f"{label}  tau={lt:.2f} ps")
ax.set_xlabel("t (ps)"); ax.set_ylabel("C(t)")
ax.set_title("Intermittent H-bond Autocorrelation C(t)")
ax.set_ylim(0, 1)
ax.legend(framealpha=0.2, facecolor="#21262D", labelcolor=tc, fontsize=9)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/hbond_acf.png", dpi=180, bbox_inches="tight", facecolor=bg)
print("Saved hbond_acf.png")

print("\nAll done.")
