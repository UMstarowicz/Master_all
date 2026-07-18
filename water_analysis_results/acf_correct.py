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
OUT    = "water_analysis_results/hbond_acf_corrected.png"

# ── KEY FIX: stride=1 → 1 ps resolution, captures fast decay ─────────
STRIDE  = 1
MAX_TAU = 300   # 300 ps is enough to see full decay

bg, axc, tc = "#0D1117", "#161B22", "#C9D1D9"

def style_ax(ax):
    ax.set_facecolor(axc)
    ax.tick_params(colors=tc, labelsize=9)
    ax.xaxis.label.set_color(tc); ax.yaxis.label.set_color(tc)
    ax.title.set_color("#E6EDF3")
    for sp in ax.spines.values(): sp.set_color("#30363D")

def compute_acf(hba, u, max_tau, stride):
    hb = hba.results.hbonds
    if len(hb) == 0:
        return np.array([]), np.array([])
    frames = hb[:,0].astype(int)
    donors = hb[:,1].astype(int)
    acc    = hb[:,3].astype(int)
    uniq   = np.unique(frames)
    n      = len(uniq)
    f2i    = {f: i for i, f in enumerate(uniq)}
    pairs  = {}
    for fr, d, a in zip(frames, donors, acc):
        k = (int(d), int(a))
        pairs.setdefault(k, set()).add(f2i[int(fr)])
    tau_max = min(max_tau, n // 2)
    num = np.zeros(tau_max)
    den = 0.0
    for ps in pairs.values():
        h = np.zeros(n, dtype=np.float32)
        for i in ps:
            h[i] = 1.0
        if h[0] == 0:
            continue
        den += 1.0
        for tau in range(tau_max):
            num[tau] += h[tau]
    if den == 0:
        return np.array([]), np.array([])
    acf = num / den
    acf /= acf[0]
    dt   = u.trajectory.dt * stride
    return np.arange(tau_max) * dt, acf

# ── Run HBA and compute ACF for each model ────────────────────────────
results = {}
for label, paths in MODELS.items():
    print(f"\n{'='*50}\n  {label}  (stride={STRIDE})\n{'='*50}")
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
    tau_ps, acf = compute_acf(hba, u, MAX_TAU, STRIDE)
    lt = float(trapezoid(acf, tau_ps)) if len(acf) > 1 else float("nan")
    print(f"  tau_HB = {lt:.2f} ps")
    results[label] = (tau_ps, acf, lt)

# ── Plot: two panels — full decay + zoom 0-30 ps ─────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.patch.set_facecolor(bg)
for ax in axes: style_ax(ax)

for label, (tau_ps, acf, lt) in results.items():
    c = COLORS[label]
    lbl = f"{label}  τ={lt:.1f} ps"
    axes[0].plot(tau_ps, acf, color=c, lw=1.8, label=lbl)
    mask = tau_ps <= 30
    axes[1].plot(tau_ps[mask], acf[mask], color=c, lw=1.8, label=lbl)

for ax, title, xlim in zip(
        axes,
        ["Intermittent C(t) — full decay (0–300 ps)",
         "Zoom: fast librational component (0–30 ps)"],
        [(0, 300), (0, 30)]):
    ax.axhline(1/np.e, color="#484F58", lw=0.8, ls=":", label="1/e ≈ 0.368")
    ax.set_xlabel("t (ps)"); ax.set_ylabel("C(t)")
    ax.set_title(title, fontsize=10)
    ax.set_ylim(0, 1); ax.set_xlim(xlim)
    ax.legend(framealpha=0.2, facecolor="#21262D", labelcolor=tc, fontsize=8)

fig.suptitle("H-bond Intermittent ACF  [stride = 1 ps — correct resolution]",
             color="#E6EDF3", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(OUT, dpi=180, bbox_inches="tight", facecolor=bg)
print(f"\nSaved → {OUT}")
