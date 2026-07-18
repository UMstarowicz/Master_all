import os, warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid as trapz
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.lib.distances import distance_array, calc_angles, calc_bonds
warnings.filterwarnings("ignore")

MODELS = {
    "FBA/eps":     {"tpr": "fbaeps/npt_correct/npt.tpr",   "xtc": "fbaeps/npt_correct/npt.xtc"},
    "SPC/Fw":      {"tpr": "spcfw/npt/npt_fixed.tpr",      "xtc": "spcfw/npt/npt.xtc"},
    "TIP4P/2005f": {"tpr": "tip4p2005f/npt/npt.tpr",       "xtc": "tip4p2005f/npt/npt.xtc"},
}
COLORS  = {"FBA/eps": "#2196F3", "SPC/Fw": "#E91E63", "TIP4P/2005f": "#4CAF50"}
OUT_DIR = "water_analysis_results"
os.makedirs(OUT_DIR, exist_ok=True)
HBOND_DIST  = 3.5
HBOND_ANGLE = 150.0
STRIDE_HBOND = 10
STRIDE_GEOM  = 5

bg  = "#0D1117"
axc = "#161B22"
tc  = "#C9D1D9"

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
    oh_lengths, hoh_angles = [], []
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

def get_hbond_geometry(u):
    oo_dists, oho_angles = [], []
    for ts in u.trajectory[::STRIDE_HBOND]:
        box  = ts.dimensions
        box3 = box[:3]
        for mol_i in u.residues:
            ow_i = mol_i.atoms.select_atoms("name OW")
            hw_i = mol_i.atoms.select_atoms("name HW*")
            if len(ow_i) == 0 or len(hw_i) == 0:
                continue
            o_donor = ow_i.positions[0]
            for h_pos in hw_i.positions:
                oh_vec = h_pos - o_donor
                oh_vec -= np.round(oh_vec / box3) * box3
                oh_len  = np.linalg.norm(oh_vec)
                if oh_len < 0.5:
                    continue
                for mol_j in u.residues:
                    if mol_i == mol_j:
                        continue
                    ow_j = mol_j.atoms.select_atoms("name OW")
                    if len(ow_j) == 0:
                        continue
                    o_acc   = ow_j.positions[0]
                    ha_vec  = o_acc - h_pos
                    ha_vec -= np.round(ha_vec / box3) * box3
                    ha_len  = np.linalg.norm(ha_vec)
                    if ha_len > 2.5:
                        continue
                    oo_vec  = o_acc - o_donor
                    oo_vec -= np.round(oo_vec / box3) * box3
                    oo_len  = np.linalg.norm(oo_vec)
                    if oo_len > HBOND_DIST:
                        continue
                    cos_ang = np.dot(oh_vec/oh_len, ha_vec/ha_len)
                    ang_deg = np.degrees(np.arccos(np.clip(cos_ang, -1, 1)))
                    if ang_deg < HBOND_ANGLE:
                        continue
                    oo_dists.append(oo_len)
                    oho_angles.append(ang_deg)
    u.trajectory.rewind()
    return np.array(oo_dists), np.array(oho_angles)

def hbond_intermittent_acf(u):
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
    tau, acf = hba.lifetime(intermittent=True, window_step=1, tau_max=200)
    u.trajectory.rewind()
    dt = u.trajectory.dt
    return np.array(tau) * STRIDE_HBOND * dt, np.array(acf)

# ── main loop ──────────────────────────────────────────────────────────
results = {}
for label, paths in MODELS.items():
    print(f"\n{'='*50}\n  {label}\n{'='*50}")
    u = mda.Universe(paths["tpr"], paths["xtc"])

    print("[1/4] RDF O-H and H-H ...")
    r_oh, g_oh = read_rdf(u, "name OW", "name HW*", rng=(0.5,5.0), excl=(1,1))
    r_hh, g_hh = read_rdf(u, "name HW*", "name HW*", rng=(0.5,5.0), excl=(1,1))

    print("[2/4] Intramolecular geometry ...")
    oh_len, hoh_ang = get_oh_hoh_geometry(u)

    print("[3/4] H-bond geometry ...")
    oo_d, oho_a = get_hbond_geometry(u)
    print(f"  Found {len(oo_d)} H-bond instances")

    print("[4/4] Intermittent H-bond ACF ...")
    try:
        tau_ps, acf = hbond_intermittent_acf(u)
        hb_life = trapz(acf, tau_ps) if len(acf) > 1 else float('nan')
        print(f"  tau_HB = {hb_life:.2f} ps")
    except Exception as e:
        print(f"  ACF failed: {e}")
        tau_ps, acf, hb_life = np.array([]), np.array([]), float('nan')

    results[label] = {
        "rdf_oh": (r_oh, g_oh), "rdf_hh": (r_hh, g_hh),
        "oh_len": oh_len, "hoh_ang": hoh_ang,
        "oo_dist": oo_d, "oho_angle": oho_a,
        "hb_acf": (tau_ps, acf), "hb_life": hb_life,
    }

# ── Figure 1: RDF O-H and H-H ─────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.patch.set_facecolor(bg)
for ax, key, title in zip(axes,
        ["rdf_oh","rdf_hh"],
        ["O-H Radial Distribution Function","H-H Radial Distribution Function"]):
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
    bins = np.linspace(0.85, 1.15, 120)
    h, e = np.histogram(d, bins=bins, density=True)
    axes[0].plot(0.5*(e[:-1]+e[1:]), h, color=c, lw=1.8,
                 label=f"{label}  mean={d.mean():.4f} A")
    a = res["hoh_ang"]
    bins = np.linspace(90, 130, 120)
    h, e = np.histogram(a, bins=bins, density=True)
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
fig.savefig(f"{OUT_DIR}/intramolecular_geometry.png", dpi=180, bbox_inches="tight", facecolor=bg)
print("Saved intramolecular_geometry.png")

# ── Figure 3: H-bond geometry ─────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.patch.set_facecolor(bg)
for ax in axes: style_ax(ax)
for label, res in results.items():
    c = COLORS[label]
    oo = res["oo_dist"]; ag = res["oho_angle"]
    if len(oo) == 0: continue
    h, e = np.histogram(oo, bins=np.linspace(2.4, 3.6, 100), density=True)
    axes[0].plot(0.5*(e[:-1]+e[1:]), h, color=c, lw=1.8,
                 label=f"{label}  mean={oo.mean():.3f} A")
    h, e = np.histogram(ag, bins=np.linspace(140, 180, 100), density=True)
    axes[1].plot(0.5*(e[:-1]+e[1:]), h, color=c, lw=1.8,
                 label=f"{label}  mean={ag.mean():.1f} deg")
axes[0].set_xlabel("O...O distance (A)"); axes[0].set_ylabel("Probability density")
axes[0].set_title("H-bond O...O Distance")
axes[0].legend(framealpha=0.2, facecolor="#21262D", labelcolor=tc, fontsize=8)
axes[1].set_xlabel("O-H...O angle (deg)"); axes[1].set_ylabel("Probability density")
axes[1].set_title("H-bond O-H...O Angle")
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
