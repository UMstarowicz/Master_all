import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import warnings
warnings.filterwarnings("ignore")

MODELS = {
    "FBA/eps":     {"density": "fbaeps/npt_correct/density.xvg",
                    "temp":    "fbaeps/npt_correct/temp.xvg"},
    "SPC/Fw":      {"density": "spcfw/npt/density.xvg",
                    "temp":    "spcfw/npt/temp.xvg"},
    "TIP4P/2005f": {"density": "tip4p2005f/npt/density.xvg",
                    "temp":    "tip4p2005f/npt/temp.xvg"},
}
COLORS = {"FBA/eps": "#2196F3", "SPC/Fw": "#E91E63", "TIP4P/2005f": "#4CAF50"}
OUT    = "water_analysis_results/convergence_distributions.png"
CUTOFF_NS = 5.0

bg, axc, tc = "#0D1117", "#161B22", "#C9D1D9"

def style_ax(ax, title, xlabel):
    ax.set_facecolor(axc)
    ax.tick_params(colors=tc, labelsize=9)
    ax.xaxis.label.set_color(tc); ax.yaxis.label.set_color(tc)
    ax.title.set_color("#E6EDF3")
    for sp in ax.spines.values(): sp.set_color("#30363D")
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(xlabel); ax.set_ylabel("Probability density")
    ax.autoscale(axis='y')

def read_xvg(path, cutoff_ns=CUTOFF_NS):
    t, v = [], []
    with open(path) as f:
        for line in f:
            if line.startswith(('#', '@')): continue
            cols = line.split()
            if len(cols) >= 2:
                t.append(float(cols[0]) / 1000.0)   # ps → ns
                v.append(float(cols[1]))
    t, v = np.array(t), np.array(v)
    mask = t <= cutoff_ns
    return v[mask]

def smooth_hist(data, n_bins=150, sigma=2.5):
    h, e = np.histogram(data, bins=n_bins, density=True)
    x = 0.5 * (e[:-1] + e[1:])

    y = gaussian_filter1d(h.astype(float), sigma=sigma)
    y = y / np.trapz(y, x)

    return x, y

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.patch.set_facecolor(bg)

# ── Density panel ─────────────────────────────────────────────────────
ax = axes[0]
style_ax(ax, "Density Distribution (0–5 ns)", "Density (kg/m³)")
for label, paths in MODELS.items():
    c    = COLORS[label]
    data = read_xvg(paths["density"])
    mean = data.mean()
    std  = data.std()
    x, y = smooth_hist(data, n_bins=120)
    ax.plot(x, y, color=c, lw=2.2,
            label=f"{label}  {mean:.1f}±{std:.1f} kg/m³")
    ax.axvline(mean, color=c, lw=1.5, ls="--")

ax.axvline(997.0, color="#FFC107", lw=1.2, ls=":",
           label="Exp. 997.0 kg/m³")
ax.legend(framealpha=0.2, facecolor="#21262D", labelcolor=tc, fontsize=9)

# ── Temperature panel ─────────────────────────────────────────────────
ax = axes[1]
style_ax(ax, "Temperature Distribution (0–5 ns)", "Temperature (K)")
for label, paths in MODELS.items():
    c    = COLORS[label]
    data = read_xvg(paths["temp"])
    mean = data.mean()
    std  = data.std()
    x, y = smooth_hist(data, n_bins=120)
    ax.plot(x, y, color=c, lw=2.2,
            label=f"{label}  {mean:.2f}±{std:.2f} K")
    ax.axvline(mean, color=c, lw=1.5, ls="--")

ax.axvline(300.0, color="#FFC107", lw=1.2, ls=":",
           label="Target 300 K")
ax.legend(framealpha=0.2, facecolor="#21262D", labelcolor=tc, fontsize=9)

fig.suptitle("NPT Ensemble Distributions — FBA/eps · SPC/Fw · TIP4P/2005f",
             color="#E6EDF3", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(OUT, dpi=180, bbox_inches="tight", facecolor=bg)
print(f"Saved → {OUT}")
