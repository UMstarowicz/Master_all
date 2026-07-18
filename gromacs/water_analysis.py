"""
Water Model Analysis Suite
==========================
Computes for FBA/epsilon, SPC/Fw, TIP4P/2005f:
  1. RDF (O–O and O–H)
  2. Self-diffusion coefficient (Einstein MSD, linear regime fit)
  3. Tetrahedral order parameter (q_tet)
  4. Coordination number (integral of RDF to first minimum)
  5. Density stability over time

Requires: MDAnalysis >= 2.0, numpy, scipy, matplotlib
"""

import os
import warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats, signal
from scipy.optimize import curve_fit
import MDAnalysis as mda
from MDAnalysis.analysis.msd import EinsteinMSD
from MDAnalysis.analysis.rdf import InterRDF

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────
# CONFIGURATION  –– edit paths and parameters
# ─────────────────────────────────────────────
MODELS = {
    "FBA/eps":       {"tpr": "fbaeps/npt_correct/npt.tpr", "xtc": "fbaeps/npt_correct/npt_nojump.xtc"},
    "SPC/Fw":        {"tpr": "spcfw/npt/npt_fixed.tpr",       "xtc": "spcfw/npt/npt_nojump.xtc"},
    "TIP4P/2005f":   {"tpr": "tip4p2005f/npt/npt.tpr",        "xtc": "tip4p2005f/npt/npt_nojump.xtc"},
}

# Water molecular mass (g/mol) – used for density if not from energy file
WATER_MASS_G_MOL = 18.015

# MSD linear-fit window: use middle [lo_frac … hi_frac] of the timeseries
MSD_FIT_LO = 0.10
MSD_FIT_HI = 0.90

# RDF parameters
RDF_NBINS  = 300
RDF_RANGE  = (0.0, 8.0)   # Å

# Tetrahedral order: how many neighbours to consider
QTET_N_NEIGHBOURS = 4
# Analyse every N-th frame (1 = all) to control runtime
QTET_STRIDE = 5

OUT_DIR = "water_analysis_results"
os.makedirs(OUT_DIR, exist_ok=True)

# ─────────────────────────────────────────────
# COLOUR PALETTE
# ─────────────────────────────────────────────
COLORS = {"FBA/ε": "#2196F3", "SPC/Fw": "#E91E63", "TIP4P/2005f": "#4CAF50"}

# ─────────────────────────────────────────────
# HELPER UTILITIES
# ─────────────────────────────────────────────

def mol_mass_from_universe(u):
    """Total molecular mass in grams from universe masses (u → g)."""
    masses_amu = u.atoms.masses          # atomic mass units (Da)
    total_amu  = masses_amu.sum()
    return total_amu / 6.02214076e23     # grams


def box_volume_A3(ts):
    """Box volume in Å³ from a timestep."""
    dims = ts.dimensions                 # [a, b, c, α, β, γ]
    # For orthorhombic/cubic boxes this is exact; good approx otherwise
    return dims[0] * dims[1] * dims[2]


def build_time_axis(msd_obj, n_frames, dt_ps=None):
    """Return time array in ps matching MSD timeseries length."""
    if dt_ps is None:
        # Fall back: infer from n_frames vs total time stored in times attr
        try:
            times = msd_obj.results.times
            dt_ps = float(times[1] - times[0]) if len(times) > 1 else 1.0
        except Exception:
            dt_ps = 1.0   # last resort: 1 ps (matches nstxout-compressed=1000, dt=0.001)
    return np.arange(n_frames) * dt_ps


# ─────────────────────────────────────────────
# 1. RDF
# ─────────────────────────────────────────────

def compute_rdf(u, label):
    """O–O RDF. Returns (r, g_r) arrays."""
    ow = u.select_atoms("name OW")
    rdf = InterRDF(ow, ow,
                   nbins=RDF_NBINS,
                   range=RDF_RANGE,
                   exclusion_block=(1, 1))   # exclude self-pair
    rdf.run(verbose=True)
    r   = rdf.results.bins          # bin centres, Å
    g_r = rdf.results.rdf
    return r, g_r


# ─────────────────────────────────────────────
# 2. SELF-DIFFUSION COEFFICIENT
# ─────────────────────────────────────────────

def compute_diffusion(u, label):
    """
    EinsteinMSD → linear fit in the diffusive regime.
    Returns D in m²/s, plus (t_ps, msd_A2) for plotting.
    """
    ow  = u.select_atoms("name OW")
    msd = EinsteinMSD(ow, select="all", msd_type="xyz", fft=True)
    msd.run(verbose=True)

    ts_data = msd.results.timeseries          # shape (N,)  Å²
    n       = len(ts_data)
    t_ps    = build_time_axis(msd, n)

    # Linear fit on the middle portion (avoid ballistic & noise)
    lo = int(n * MSD_FIT_LO)
    hi = int(n * MSD_FIT_HI)
    slope, intercept, r, p, se = stats.linregress(t_ps[lo:hi], ts_data[lo:hi])

    # D = slope / 6   (3-D MSD: <r²> = 6Dt)
    # slope in Å²/ps → convert to m²/s:  1 Å²/ps = 1e-20 / 1e-12 = 1e-8 m²/s
    D_m2s = (slope / 6.0) * 1e-8
    print(f"  [{label}] D = {D_m2s:.3e} m²/s  (R²={r**2:.4f}, slope={slope:.4f} Å²/ps)")

    return D_m2s, t_ps, ts_data, slope, intercept, lo, hi


# ─────────────────────────────────────────────
# 3. TETRAHEDRAL ORDER PARAMETER
# ─────────────────────────────────────────────

def _qtet_frame(positions, box):
    """
    q_tet for one frame using distance_array (no double-counting).
    q = 1 - (3/8) Σ_{j<k} (cos θ_{jk} + 1/3)²  averaged over molecules.
    """
    from MDAnalysis.lib.distances import distance_array

    n = len(positions)
    # Full N×N distance matrix with PBC minimum image
    dmat = distance_array(positions, positions, box=box)  # shape (n, n)

    q_sum   = 0.0
    counted = 0
    box3    = box[:3]

    for i in range(n):
        # Get 4 nearest neighbours (exclude self by setting diagonal to inf)
        row = dmat[i].copy()
        row[i] = np.inf
        nn_idx = np.argpartition(row, QTET_N_NEIGHBOURS)[:QTET_N_NEIGHBOURS]
        # Verify they are within a reasonable shell (< 4 Å)
        if row[nn_idx].max() > 4.0:
            continue

        # Compute minimum-image unit vectors from i to each neighbour
        vecs = []
        for j in nn_idx:
            dr   = positions[j] - positions[i]
            dr  -= np.round(dr / box3) * box3   # minimum image
            norm = np.linalg.norm(dr)
            if norm > 0:
                vecs.append(dr / norm)

        if len(vecs) < QTET_N_NEIGHBOURS:
            continue

        # Sum over all unique pairs of neighbours
        q_i = 0.0
        for a in range(QTET_N_NEIGHBOURS):
            for b in range(a + 1, QTET_N_NEIGHBOURS):
                cos_ab = np.clip(np.dot(vecs[a], vecs[b]), -1.0, 1.0)
                q_i   += (cos_ab + 1.0 / 3.0) ** 2
        q_sum  += 1.0 - (3.0 / 8.0) * q_i
        counted += 1

    return q_sum / counted if counted > 0 else 0.0


def compute_qtet(u, label):
    """
    Tetrahedral order parameter per frame (strided).
    Returns (frame_times_ps, q_array).
    """
    ow  = u.select_atoms("name OW")
    qs  = []
    ts_list = []
    total = len(u.trajectory[::QTET_STRIDE])
    print(f"  [{label}] q_tet: processing {total} frames …")
    for k, ts in enumerate(u.trajectory[::QTET_STRIDE]):
        if k % max(1, total // 10) == 0:
            print(f"    {k}/{total}", flush=True)
        q = _qtet_frame(ow.positions.copy(), ts.dimensions)
        qs.append(q)
        ts_list.append(ts.time)
    q_arr = np.array(qs)
    print(f"  [{label}] <q_tet> = {q_arr.mean():.4f} ± {q_arr.std():.4f}")
    return np.array(ts_list), q_arr


# ─────────────────────────────────────────────
# 4. COORDINATION NUMBER
# ─────────────────────────────────────────────

def compute_coordination_number(r, g_r, rho_bulk, r_min=2.4, r_max_search_start=3.0):
    """
    Integrate g(r)*4πr² ρ up to the first minimum after r_min.
    rho_bulk: bulk O density in Å⁻³.
    """
    # Find first minimum of g_r between r_min and some upper limit
    mask = (r > r_min) & (r < 5.0)
    g_sub = g_r[mask]
    r_sub = r[mask]
    # Local minima
    minima, _ = signal.find_peaks(-g_sub)
    if len(minima) == 0:
        r_cutoff = 3.5   # fallback
    else:
        r_cutoff = r_sub[minima[0]]

    # Integrate
    mask_int = r <= r_cutoff
    integrand = g_r[mask_int] * 4.0 * np.pi * r[mask_int]**2
    cn = np.trapz(integrand, r[mask_int]) * rho_bulk
    print(f"  Coordination number (r_cut={r_cutoff:.2f} Å): {cn:.2f}")
    return cn, r_cutoff


# ─────────────────────────────────────────────
# 5. DENSITY OVER TIME
# ─────────────────────────────────────────────

def compute_density(u, label):
    """
    Density in g/cm³ from box volume at each frame.
    Uses atom masses from the TPR.
    """
    total_g = mol_mass_from_universe(u)
    densities = []
    times     = []
    for ts in u.trajectory:
        vol_A3  = box_volume_A3(ts)
        vol_cm3 = vol_A3 * 1e-24           # 1 Å³ = 1e-24 cm³
        densities.append(total_g / vol_cm3)
        times.append(ts.time)
    rho   = np.array(densities)
    t_ps  = np.array(times)
    print(f"  [{label}] <ρ> = {rho.mean():.4f} ± {rho.std():.4f} g/cm³")
    return t_ps, rho


# ─────────────────────────────────────────────
# PLOTTING
# ─────────────────────────────────────────────

def plot_all(results):
    """Master figure: 5 panels, one curve per model."""
    fig = plt.figure(figsize=(18, 14))
    fig.patch.set_facecolor("#0D1117")

    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.35)
    axes = [
        fig.add_subplot(gs[0, 0]),   # RDF
        fig.add_subplot(gs[0, 1]),   # MSD
        fig.add_subplot(gs[1, 0]),   # q_tet
        fig.add_subplot(gs[1, 1]),   # Coord. number (bar)
        fig.add_subplot(gs[2, :]),   # Density
    ]

    panel_style = dict(facecolor="#161B22")
    for ax in axes:
        ax.set_facecolor("#161B22")
        ax.tick_params(colors="#C9D1D9", labelsize=9)
        ax.xaxis.label.set_color("#C9D1D9")
        ax.yaxis.label.set_color("#C9D1D9")
        ax.title.set_color("#E6EDF3")
        for spine in ax.spines.values():
            spine.set_color("#30363D")

    # ── Panel 1: RDF ──────────────────────────────────────────
    ax = axes[0]
    for label, res in results.items():
        r, g_r = res["rdf"]
        ax.plot(r, g_r, color=COLORS[label], lw=1.8, label=label)
    ax.axhline(1.0, color="#484F58", lw=0.8, ls="--")
    ax.set_xlabel("r (Å)")
    ax.set_ylabel("g(r)")
    ax.set_title("O–O Radial Distribution Function")
    ax.set_xlim(2.0, 8.0)
    ax.legend(framealpha=0.2, facecolor="#21262D", labelcolor="#C9D1D9", fontsize=8)

    # ── Panel 2: MSD ──────────────────────────────────────────
    ax = axes[1]
    for label, res in results.items():
        t, msd_data, slope, intercept, lo, hi = res["msd_plot"]
        ax.plot(t, msd_data, color=COLORS[label], lw=1.4, alpha=0.7, label=label)
        # Show fit
        t_fit = t[lo:hi]
        ax.plot(t_fit, slope * t_fit + intercept,
                color=COLORS[label], lw=2.0, ls="--", alpha=1.0)
    ax.set_xlabel("t (ps)")
    ax.set_ylabel("MSD (Å²)")
    ax.set_title("Mean Squared Displacement  (dashed = linear fit)")
    ax.legend(framealpha=0.2, facecolor="#21262D", labelcolor="#C9D1D9", fontsize=8)

    # ── Panel 3: q_tet over time ──────────────────────────────
    ax = axes[2]
    for label, res in results.items():
        t_q, q_arr = res["qtet"]
        ax.plot(t_q / 1000, q_arr, color=COLORS[label], lw=1.2, alpha=0.6)
        ax.axhline(q_arr.mean(), color=COLORS[label], lw=1.8, ls="-",
                   label=f"{label}  ⟨q⟩={q_arr.mean():.3f}")
    ax.set_xlabel("t (ns)")
    ax.set_ylabel("q_tet")
    ax.set_title("Tetrahedral Order Parameter")
    ax.set_ylim(0, 1)
    ax.legend(framealpha=0.2, facecolor="#21262D", labelcolor="#C9D1D9", fontsize=8)

    # ── Panel 4: Coordination number bar ─────────────────────
    ax = axes[3]
    labels_bar = list(results.keys())
    cn_vals    = [results[l]["cn"] for l in labels_bar]
    colors_bar = [COLORS[l] for l in labels_bar]
    bars = ax.bar(labels_bar, cn_vals, color=colors_bar, alpha=0.85, width=0.5,
                  edgecolor="#30363D", linewidth=1.2)
    ax.axhline(4.0, color="#FFC107", lw=1.2, ls="--", label="Ideal = 4.0")
    for bar, val in zip(bars, cn_vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.05,
                f"{val:.2f}", ha="center", va="bottom",
                color="#E6EDF3", fontsize=10, fontweight="bold")
    ax.set_ylabel("Coordination number")
    ax.set_title("First-shell Coordination Number")
    ax.legend(framealpha=0.2, facecolor="#21262D", labelcolor="#C9D1D9", fontsize=8)
    ax.set_ylim(0, max(cn_vals) * 1.25)

    # ── Panel 5: Density over time ────────────────────────────
    ax = axes[4]
    for label, res in results.items():
        t_rho, rho = res["density"]
        ax.plot(t_rho / 1000, rho, color=COLORS[label], lw=1.2, alpha=0.55)
        # Rolling mean
        win = max(1, len(rho) // 50)
        rho_smooth = np.convolve(rho, np.ones(win) / win, mode="valid")
        t_smooth   = t_rho[win - 1:] / 1000
        ax.plot(t_smooth, rho_smooth, color=COLORS[label], lw=2.2,
                label=f"{label}  ⟨ρ⟩={rho.mean():.4f} g/cm³")
    ax.axhline(0.9970, color="#FFC107", lw=1.0, ls=":", label="Exp 298 K (0.9970)")
    ax.set_xlabel("t (ns)")
    ax.set_ylabel("ρ (g/cm³)")
    ax.set_title("Density Stability Over Time  (thin=raw, thick=rolling mean)")
    ax.legend(framealpha=0.2, facecolor="#21262D", labelcolor="#C9D1D9",
              fontsize=8, ncol=2)

    # Global title
    fig.suptitle("Water Model Comparison: FBA/ε · SPC/Fw · TIP4P/2005f",
                 color="#E6EDF3", fontsize=14, fontweight="bold", y=0.98)

    out = os.path.join(OUT_DIR, "water_model_analysis.png")
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(f"\nFigure saved → {out}")
    return out


def save_summary_csv(results):
    """Write a tidy CSV with scalar observables."""
    rows = ["model,D_m2s,q_tet_mean,q_tet_std,coord_number,rho_mean_gcm3,rho_std_gcm3"]
    for label, res in results.items():
        D     = res["D"]
        q_m   = res["qtet"][1].mean()
        q_s   = res["qtet"][1].std()
        cn    = res["cn"]
        rho_m = res["density"][1].mean()
        rho_s = res["density"][1].std()
        rows.append(f"{label},{D:.4e},{q_m:.4f},{q_s:.4f},{cn:.4f},{rho_m:.6f},{rho_s:.6f}")
    out = os.path.join(OUT_DIR, "summary.csv")
    with open(out, "w") as f:
        f.write("\n".join(rows))
    print(f"Summary CSV → {out}")


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    results = {}

    for label, paths in MODELS.items():
        tpr, xtc = paths["tpr"], paths["xtc"]
        if not (os.path.exists(tpr) and os.path.exists(xtc)):
            print(f"\n[SKIP] {label}: files not found ({tpr}, {xtc})")
            continue

        print(f"\n{'='*60}")
        print(f"  Analysing {label}")
        print(f"{'='*60}")

        u = mda.Universe(tpr, xtc)
        n_ow    = u.select_atoms("name OW").n_atoms
        n_water = n_ow                        # one OW per molecule
        print(f"  Frames: {u.trajectory.n_frames}  |  OW atoms: {n_ow}")

        # 1. RDF
        print("\n[1/5] RDF …")
        r, g_r = compute_rdf(u, label)
        u.trajectory.rewind()

        # Bulk O number density (atoms/Å³) from average volume
        vol_list = []
        for ts in u.trajectory[::10]:
            vol_list.append(box_volume_A3(ts))
        rho_ow = n_ow / np.mean(vol_list)
        u.trajectory.rewind()

        # 2. Diffusion
        print("\n[2/5] Self-diffusion …")
        D, t_msd, msd_data, slope, intercept, lo, hi = compute_diffusion(u, label)
        u.trajectory.rewind()

        # 3. q_tet
        print("\n[3/5] Tetrahedral order …")
        t_q, q_arr = compute_qtet(u, label)
        u.trajectory.rewind()

        # 4. Coordination number
        print("\n[4/5] Coordination number …")
        cn, r_cut = compute_coordination_number(r, g_r, rho_bulk=rho_ow)

        # 5. Density
        print("\n[5/5] Density over time …")
        t_rho, rho = compute_density(u, label)
        u.trajectory.rewind()

        results[label] = {
            "rdf":      (r, g_r),
            "D":        D,
            "msd_plot": (t_msd, msd_data, slope, intercept, lo, hi),
            "qtet":     (t_q, q_arr),
            "cn":       cn,
            "density":  (t_rho, rho),
        }

    if not results:
        print("\nNo trajectories found. Check MODELS paths at the top of the script.")
        return

    # ── Summary table ──────────────────────────────────────────
    print(f"\n{'='*60}")
    print("  SUMMARY")
    print(f"{'='*60}")
    print(f"{'Model':<20} {'D (m²/s)':>14} {'q_tet':>8} {'CN':>6} {'ρ (g/cm³)':>12}")
    print("-" * 65)
    for label, res in results.items():
        q   = res["qtet"][1].mean()
        rho = res["density"][1].mean()
        print(f"{label:<20} {res['D']:>14.3e} {q:>8.4f} {res['cn']:>6.2f} {rho:>12.4f}")

    # ── Plots & CSV ─────────────────────────────────────────────
    plot_all(results)
    save_summary_csv(results)


if __name__ == "__main__":
    main()
