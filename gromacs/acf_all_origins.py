import warnings, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
warnings.filterwarnings("ignore")

MODELS = {
    "FBA/eps":     {"tpr": "fbaeps/npt_correct/npt.tpr",  "xtc": "fbaeps/npt_correct/npt.xtc"},
    "SPC/Fw":      {"tpr": "spcfw/npt/npt_fixed.tpr",     "xtc": "spcfw/npt/npt.xtc"},
    "TIP4P/2005f": {"tpr": "tip4p2005f/npt/npt.tpr",      "xtc": "tip4p2005f/npt/npt.xtc"},
}
COLORS = {"FBA/eps": "#2196F3", "SPC/Fw": "#E91E63", "TIP4P/2005f": "#4CAF50"}
OUT    = "water_analysis_results/hbond_acf_final.png"
STRIDE  = 1       # 1 ps resolution
MAX_TAU = 100     # 100 ps window

bg, axc, tc = "#0D1117", "#161B22", "#C9D1D9"

def style_ax(ax):
    ax.set_facecolor(axc)
    ax.tick_params(colors=tc, labelsize=9)
    ax.xaxis.label.set_color(tc); ax.yaxis.label.set_color(tc)
    ax.title.set_color("#E6EDF3")
    for sp in ax.spines.values(): sp.set_color("#30363D")

def compute_acf_all_origins(hba, u, max_tau, stride):
    """
    Correct intermittent ACF averaged over ALL time origins.
    C(tau) = sum_t [ h(t)*h(t+tau) ] / sum_t [ h(t) ]
    This is the proper Luzar-Chandler definition.
    """
    hb = hba.results.hbonds
    if len(hb) == 0:
        return np.array([]), np.array([])

    frames   = hb[:, 0].astype(int)
    donors   = hb[:, 1].astype(int)
    acceptors= hb[:, 3].astype(int)

    uniq_frames = np.unique(frames)
    n_frames    = len(uniq_frames)
    f2i         = {f: i for i, f in enumerate(uniq_frames)}

    # Build binary presence matrix: pairs × frames
    # Use dict of arrays for memory efficiency
    pair_h = {}
    for fr, d, a in zip(frames, donors, acceptors):
        key = (int(d), int(a))
        if key not in pair_h:
            pair_h[key] = np.zeros(n_frames, dtype=np.float32)
        pair_h[key][f2i[int(fr)]] = 1.0

    tau_max = min(max_tau, n_frames // 4)
    num = np.zeros(tau_max)
    den = np.zeros(tau_max)

    print(f"  Computing ACF over {len(pair_h)} pairs, {n_frames} frames, "
          f"max_tau={tau_max}...")

    for key, h in pair_h.items():
        # Average over ALL time origins t0
        for tau in range(tau_max):
            valid = n_frames - tau
            if valid <= 0:
                break
            num[tau] += np.dot(h[:valid], h[tau:tau+valid])
            den[tau] += h[:valid].sum()

    # Avoid division by zero
    mask = den > 0
    acf  = np.zeros(tau_max)
    acf[mask] = num[mask] / den[mask]

    # Normalise so C(0) = 1
    if acf[0] > 0:
        acf /= acf[0]

    dt = u.trajectory.dt * stride
    return np.arange(tau_max) * dt, acf

# ── Main loop ─────────────────────────────────────────────────────────
results = {}
for label, paths in MODELS.items():
    print(f"\n{'='*50}\n  {label}\n{'='*50}")
    u = mda.Universe(paths["tpr"], paths["xtc"])
    hba = HydrogenBondAnalysis(
        universe=u,
        donors_sel="name OW",
        hydrogens_sel="name HW1 HW2",
        acceptors_sel="name OW",
        d_h_cutoff=1.2,
        d_a_cutoff=3.5,
        d_h_a_angle_cutoff=150.0,
        update_selections=False,
    )
    hba.run(step=STRIDE, verbose=True)
    u.trajectory.rewind()

    tau_ps, acf = compute_acf_all_origins(hba, u, MAX_TAU, STRIDE)
    lt = float(trapezoid(acf, tau_ps)) if len(acf) > 1 else float("nan")
    print(f"  tau_HB (all origins) = {lt:.2f} ps")
    results[label] = (tau_ps, acf, lt)

# ── Plot ──────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.patch.set_facecolor(bg)
for ax in axes: style_ax(ax)

for label, (tau_ps, acf, lt) in results.items():
    c = COLORS[label]
    lbl = f"{label}  τ={lt:.1f} ps"
    axes[0].plot(tau_ps, acf, color=c, lw=1.8, label=lbl)
    mask = tau_ps <= 20
    axes[1].plot(tau_ps[mask], acf[mask], color=c, lw=1.8, label=lbl)

for ax, title, xlim in zip(axes,
        ["Intermittent C(t) — full decay",
         "Zoom: 0–20 ps (fast component)"],
        [(0, 100), (0, 20)]):
    ax.axhline(1/np.e, color="#FFC107", lw=1.0, ls=":",
               label="1/e ≈ 0.368")
    ax.axhspan(0, 2, alpha=0.08, color="#4CAF50",
               label="Exp. τ range (1–2 ps)")
    ax.set_xlabel("t (ps)"); ax.set_ylabel("C(t)")
    ax.set_title(title, fontsize=10)
    ax.set_ylim(0, 1); ax.set_xlim(xlim)
    ax.legend(framealpha=0.2, facecolor="#21262D",
              labelcolor=tc, fontsize=8)

fig.suptitle(
    "Intermittent H-bond ACF  [all time origins, stride=1 ps]\n"
    "Note: residual overestimate vs exp (~1–2 ps) due to 1 ps frame spacing",
    color="#E6EDF3", fontsize=11, fontweight="bold")
fig.tight_layout()
fig.savefig(OUT, dpi=180, bbox_inches="tight", facecolor=bg)
print(f"\nSaved → {OUT}")
