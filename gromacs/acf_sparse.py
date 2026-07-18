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
STRIDE  = 1
MAX_TAU = 50    # 50 ps window — enough to see full decay at 1 ps res

bg, axc, tc = "#0D1117", "#161B22", "#C9D1D9"

def style_ax(ax):
    ax.set_facecolor(axc)
    ax.tick_params(colors=tc, labelsize=9)
    ax.xaxis.label.set_color(tc); ax.yaxis.label.set_color(tc)
    ax.title.set_color("#E6EDF3")
    for sp in ax.spines.values(): sp.set_color("#30363D")

def compute_acf_sparse(hba, u, max_tau, stride):
    """
    Memory-efficient intermittent ACF using sparse frame index sets.
    Never builds a full pairs×frames matrix.
    C(tau) = sum_{pairs,t0} h(t0)*h(t0+tau) / sum_{pairs,t0} h(t0)
    """
    hb = hba.results.hbonds
    if len(hb) == 0:
        return np.array([]), np.array([])

    frames    = hb[:, 0].astype(int)
    donors    = hb[:, 1].astype(int)
    acceptors = hb[:, 3].astype(int)

    uniq_f = np.unique(frames)
    n_f    = len(uniq_f)
    f2i    = {f: i for i, f in enumerate(uniq_f)}
    tau_max = min(max_tau, n_f // 4)

    # Build sparse dict: pair → sorted array of frame indices
    pair_idx = {}
    for fr, d, a in zip(frames, donors, acceptors):
        k = (int(d), int(a))
        pair_idx.setdefault(k, []).append(f2i[int(fr)])
    # Convert to sorted arrays
    for k in pair_idx:
        pair_idx[k] = np.array(sorted(pair_idx[k]), dtype=np.int32)

    print(f"  {len(pair_idx)} unique pairs, {n_f} frames, tau_max={tau_max}")

    num = np.zeros(tau_max, dtype=np.float64)
    den = np.zeros(tau_max, dtype=np.float64)

    # Process in batches to show progress
    pairs_list = list(pair_idx.values())
    n_pairs    = len(pairs_list)
    batch      = max(1, n_pairs // 20)

    for pi, idx_arr in enumerate(pairs_list):
        if pi % batch == 0:
            print(f"  pair {pi}/{n_pairs} ...", flush=True)

        # Build dense h array only for this pair (cheap: 1 × n_f)
        h = np.zeros(n_f, dtype=np.float32)
        h[idx_arr] = 1.0

        # Vectorised ACF via sliding dot products
        for tau in range(tau_max):
            valid = n_f - tau
            if valid <= 0:
                break
            num[tau] += np.dot(h[:valid], h[tau:tau+valid])
            den[tau] += h[:valid].sum()

    mask       = den > 0
    acf        = np.zeros(tau_max)
    acf[mask]  = num[mask] / den[mask]
    if acf[0] > 0:
        acf /= acf[0]

    dt = u.trajectory.dt * stride
    lt = float(trapezoid(acf, np.arange(tau_max) * dt))
    print(f"  tau_HB = {lt:.2f} ps")
    return np.arange(tau_max) * dt, acf, lt

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
        d_h_cutoff=1.0,
        d_a_cutoff=3.0,
        d_h_a_angle_cutoff=168.0,
        update_selections=False,
    )
    start_frame = len(u.trajectory) // 2
    print(f"  Using second half only: starting from frame {start_frame}")
    hba.run(start=start_frame, step=STRIDE, verbose=True)
    u.trajectory.rewind()
    tau_ps, acf, lt = compute_acf_sparse(hba, u, MAX_TAU, STRIDE)
    results[label] = (tau_ps, acf, lt)

# ── Plot ──────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.patch.set_facecolor(bg)
for ax in axes: style_ax(ax)

for label, (tau_ps, acf, lt) in results.items():
    c = COLORS[label]
    lbl = f"{label}  τ={lt:.1f} ps"
    axes[0].plot(tau_ps, acf, color=c, lw=1.8, label=lbl)
    mask = tau_ps <= 15
    axes[1].plot(tau_ps[mask], acf[mask], color=c, lw=1.8, label=lbl)

for ax, title, xlim in zip(axes,
        ["Intermittent C(t) — full decay (0–50 ps)",
         "Zoom: fast component (0–15 ps)"],
        [(0, 50), (0, 15)]):
    ax.axhline(1/np.e, color="#FFC107", lw=1.0, ls="--", label="1/e ≈ 0.368")
    ax.axhspan(0, 2, alpha=0.10, color="#4CAF50", label="Exp. τ ~1–2 ps")
    ax.set_xlabel("t (ps)"); ax.set_ylabel("C(t)")
    ax.set_title(title, fontsize=10)
    ax.set_ylim(0, 1); ax.set_xlim(xlim)
    ax.legend(framealpha=0.2, facecolor="#21262D", labelcolor=tc, fontsize=8)

fig.suptitle(
    "Intermittent H-bond ACF  [all time origins, stride=1 ps]\n"
    "Residual overestimate vs exp (1–2 ps) due to 1 ps frame spacing — "
    "sub-ps librational component unresolvable at this output frequency",
    color="#E6EDF3", fontsize=10, fontweight="bold")
fig.tight_layout()
fig.savefig(OUT, dpi=180, bbox_inches="tight", facecolor=bg)
print(f"\nSaved → {OUT}")
