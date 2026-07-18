"""
analyse_trajectory.py
=====================
Full structural and dynamical analysis of GROMACS NPT water trajectories.

Computes:
  - RDF:     O-O, O-H, H-H  (radial distribution functions)
  - CN:      Coordination number (integral of O-O RDF up to first minimum)
  - MSD:     Mean squared displacement → self-diffusion coefficient D
  - Density: From energy file
  - Pressure: From energy file
  - Custom:  Tetrahedral order parameter Q, H-bond statistics

Usage
-----
    # Analyse one model:
    python analyse_trajectory.py --model spcfw \
                                  --traj npt/npt.xtc \
                                  --tpr  npt/npt.tpr \
                                  --edr  npt/npt.edr \
                                  --sites 3

    # Analyse all three and compare:
    python analyse_trajectory.py --compare \
        --models spcfw fbaeps tip4p2005f \
        --dirs   retis_spcfw/npt retis_fbaeps/npt retis_tip4p2005f/npt \
        --sites  3 3 4

Dependencies
------------
    pip install MDAnalysis numpy scipy matplotlib
"""

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import argrelmin
from scipy.integrate import cumulative_trapezoid

# MDAnalysis for trajectory reading
try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rdf as mda_rdf
    from MDAnalysis.analysis.msd import EinsteinMSD
    HAS_MDA = True
except ImportError:
    HAS_MDA = False
    warnings.warn("MDAnalysis not found. Install with: pip install MDAnalysis")

KB  = 1.380649e-23     # J/K
NA  = 6.02214076e23
PS  = 1e-12            # seconds
NM2 = 1e-18            # nm² → m²
COLORS = {"spcfw": "#0072B2", "fbaeps": "#D55E00", "tip4p2005f": "#009E73"}
FALLBACK_COLORS = ["#0072B2", "#D55E00", "#009E73"]


# =========================================================================== #
#  RDF + Coordination Number                                                   #
# =========================================================================== #

def compute_rdf(universe, sel1, sel2, label, n_bins=300, range=(0.1, 0.8),
                exclusion_block=None):
    """
    Compute RDF between two atom selections using MDAnalysis.

    Parameters
    ----------
    universe         : MDAnalysis Universe
    sel1, sel2       : atom selection strings (e.g. 'name OW', 'name HW*')
    label            : string label for output
    n_bins           : number of bins
    range            : (r_min, r_max) in nm
    exclusion_block  : tuple (n, m) to exclude intramolecular pairs, e.g. (1,1)

    Returns
    -------
    r      : bin centres (nm)
    g_r    : g(r) values
    """
    if not HAS_MDA:
        raise ImportError("MDAnalysis required for RDF computation.")

    ag1 = universe.select_atoms(sel1)
    ag2 = universe.select_atoms(sel2)

    # Convert nm range to Angstrom (MDAnalysis uses Å internally)
    rng_ang = (range[0] * 10, range[1] * 10)
    bins_ang = n_bins

    rdf_obj = mda_rdf.InterRDF(
        ag1, ag2,
        nbins=bins_ang,
        range=rng_ang,
        exclusion_block=exclusion_block,
        verbose=False,
    )
    rdf_obj.run()

    r   = rdf_obj.results.bins / 10.0   # Å → nm
    g_r = rdf_obj.results.rdf
    return r, g_r


def coordination_number(r, g_r, rho, r_cut=None):
    """
    Integrate RDF to get coordination number up to first minimum.

    CN = 4π ρ ∫₀^r_cut g(r) r² dr

    Parameters
    ----------
    r      : array of distances (nm)
    g_r    : g(r)
    rho    : number density (nm^-3)
    r_cut  : integration cutoff (nm); if None, use first minimum of g(r)

    Returns
    -------
    cn      : coordination number (float)
    r_cut   : cutoff used (nm)
    """
    if r_cut is None:
        # Find first minimum after main peak
        minima = argrelmin(g_r, order=5)[0]
        # Main OO peak is around 0.28 nm; first min after that
        main_peak_idx = np.argmax(g_r[r > 0.2])
        candidates = minima[minima > (np.where(r > 0.2)[0][0] + main_peak_idx)]
        r_cut = r[candidates[0]] if len(candidates) else 0.35

    mask = r <= r_cut
    integrand = 4 * np.pi * rho * g_r[mask] * r[mask]**2
    cn = np.trapz(integrand, r[mask])
    return float(cn), float(r_cut)


# =========================================================================== #
#  MSD and Diffusion Coefficient                                               #
# =========================================================================== #

def compute_msd(universe, sel='name OW', n_frames_skip=0):
    """
    Compute MSD of oxygen atoms → self-diffusion coefficient D.

    Uses Einstein relation: D = MSD(t) / (6t)  in the linear regime.

    Returns
    -------
    times  : array (ps)
    msd    : array (nm²)
    D      : diffusion coefficient (m²/s)
    """
    if not HAS_MDA:
        raise ImportError("MDAnalysis required for MSD computation.")

    ag = universe.select_atoms(sel)
    msd_obj = EinsteinMSD(ag, select='all', msd_type='xyz', fft=True)
    msd_obj.run()

    times = msd_obj.results.timeseries[:, 0]   # ps
    msd   = msd_obj.results.timeseries[:, 1]   # Å²
    msd_nm2 = msd / 100.0                       # Å² → nm²

    # Linear fit over 20–80% of the trajectory (avoid ballistic + noise regimes)
    n = len(times)
    lo, hi = int(0.2 * n), int(0.8 * n)
    if hi - lo < 10:
        lo, hi = 1, n - 1

    t_fit   = times[lo:hi] * PS          # ps → s
    msd_fit = msd_nm2[lo:hi] * NM2        # nm² → m²

    coeffs  = np.polyfit(t_fit, msd_fit, 1)
    D       = coeffs[0] / 6.0            # m²/s  (3D: MSD = 6Dt)

    return times, msd_nm2, D


# =========================================================================== #
#  Tetrahedral Order Parameter                                                 #
# =========================================================================== #

def tetrahedral_order(universe, sel='name OW', n_neighbours=4, n_frames=100):
    """
    Compute tetrahedral order parameter Q for each oxygen.

    Q = 1 - (3/8) Σ_{j<k} (cos θ_jik + 1/3)²

    where j, k are the 4 nearest oxygen neighbours of i.
    Q = 1 for perfect tetrahedron, Q = 0 for ideal gas.

    Returns
    -------
    Q_mean : float
    Q_std  : float
    """
    from MDAnalysis.lib.distances import distance_array

    oxygens = universe.select_atoms(sel)
    Q_all   = []

    step = max(1, len(universe.trajectory) // n_frames)
    for ts in universe.trajectory[::step]:
        box   = ts.dimensions
        pos   = oxygens.positions                   # Å

        dist_matrix = distance_array(pos, pos, box=box)
        np.fill_diagonal(dist_matrix, np.inf)

        for i in range(len(oxygens)):
            nn_idx = np.argsort(dist_matrix[i])[:n_neighbours]
            nn_pos = pos[nn_idx] - pos[i]

            # Apply minimum image convention
            for dim in range(3):
                nn_pos[:, dim] -= box[dim] * np.round(nn_pos[:, dim] / box[dim])

            # Normalise
            norms = np.linalg.norm(nn_pos, axis=1, keepdims=True)
            norms[norms == 0] = 1.0
            nn_pos = nn_pos / norms

            Q_i = 0.0
            for j in range(n_neighbours):
                for k in range(j+1, n_neighbours):
                    cos_theta = np.dot(nn_pos[j], nn_pos[k])
                    Q_i += (cos_theta + 1.0/3.0)**2

            Q_all.append(1.0 - (3.0/8.0) * Q_i)

    return float(np.mean(Q_all)), float(np.std(Q_all))


# =========================================================================== #
#  Hydrogen bond statistics                                                    #
# =========================================================================== #

def hbond_stats(universe, donor_sel='name OW', acceptor_sel='name OW',
                h_sel='name HW* or name HW1 or name HW2',
                d_a_cutoff=3.5, angle_cutoff=30.0, n_frames=200):
    """
    Simple geometric H-bond counter.

    Criteria (standard):
      - D–A distance < d_a_cutoff (Å, default 3.5 Å)
      - D–H···A angle > (180 - angle_cutoff)° (default > 150°)

    Returns
    -------
    n_hbonds_mean : float — mean H-bonds per molecule
    n_hbonds_std  : float
    """
    from MDAnalysis.lib.distances import calc_angles, distance_array

    donors    = universe.select_atoms(donor_sel)
    acceptors = universe.select_atoms(acceptor_sel)
    hydrogens = universe.select_atoms(h_sel)

    hb_per_frame = []
    step = max(1, len(universe.trajectory) // n_frames)

    for ts in universe.trajectory[::step]:
        box  = ts.dimensions
        d_pos = donors.positions
        a_pos = acceptors.positions
        h_pos = hydrogens.positions

        # D–A distance matrix
        da_dist = distance_array(d_pos, a_pos, box=box)
        n_hb = 0
        for di, d in enumerate(donors):
            for ai, a in enumerate(acceptors):
                if di == ai:
                    continue
                if da_dist[di, ai] >= d_a_cutoff:
                    continue
                # Find H atoms belonging to this donor molecule
                mol_h = universe.select_atoms(
                    f"resid {d.resid} and (name HW1 or name HW2 or name HW*)"
                )
                for h in mol_h:
                    angle = calc_angles(
                        d.position, h.position, a.position, box=box
                    )
                    angle_deg = np.degrees(angle)
                    if angle_deg > (180.0 - angle_cutoff):
                        n_hb += 1
        n_molecules = len(donors)
        hb_per_frame.append(n_hb / max(n_molecules, 1))

    return float(np.mean(hb_per_frame)), float(np.std(hb_per_frame))


# =========================================================================== #
#  Energy file reader (density, pressure, temperature)                        #
# =========================================================================== #

def read_xvg(filepath):
    """Read a GROMACS .xvg file, return (times, values)."""
    times, values = [], []
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(('#', '@')) or not line:
                continue
            parts = line.split()
            try:
                times.append(float(parts[0]))
                values.append([float(x) for x in parts[1:]])
            except ValueError:
                continue
    return np.array(times), np.array(values)


def extract_energy_properties(edr_path, out_dir, model):
    """
    Use gmx energy to extract density, temperature, pressure from .edr file.
    Returns dict of {property: (mean, std)}.
    """
    import subprocess

    props = {
        "Density"    : "22",
        "Temperature": "15",
        "Pressure"   : "16",
        "Potential"  : "10",
    }
    results = {}
    for name, idx in props.items():
        xvg = Path(out_dir) / f"{model}_{name.lower()}.xvg"
        try:
            proc = subprocess.run(
                ["gmx", "energy", "-f", str(edr_path), "-o", str(xvg)],
                input=f"{idx}\n0\n", capture_output=True, text=True, timeout=60
            )
            if xvg.exists():
                t, v = read_xvg(xvg)
                if len(v):
                    # Skip first 20% as additional equilibration
                    lo = len(v) // 5
                    vals = v[lo:, 0]
                    results[name] = (float(np.mean(vals)), float(np.std(vals)))
        except Exception as e:
            warnings.warn(f"Could not extract {name}: {e}")
            results[name] = (np.nan, np.nan)
    return results


# =========================================================================== #
#  Main analysis routine                                                       #
# =========================================================================== #

def analyse_model(model, traj, tpr, edr, sites, out_dir="./analysis_results"):
    """
    Run full analysis for one model.

    Parameters
    ----------
    model   : str   — label (e.g. 'spcfw')
    traj    : str   — path to .xtc trajectory
    tpr     : str   — path to .tpr run input
    edr     : str   — path to .edr energy file
    sites   : int   — atoms per molecule (3 for SPC/Fw, FBA/eps; 4 for TIP4P/2005f)
    out_dir : str   — output directory
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  Analysing: {model}  ({sites}-site)")
    print(f"  Traj : {traj}")
    print(f"{'='*60}")

    if not HAS_MDA:
        print("ERROR: MDAnalysis not installed. Run: pip install MDAnalysis")
        return {}

    # Load trajectory
    u = mda.Universe(str(tpr), str(traj))
    n_atoms  = len(u.atoms)
    n_mol    = n_atoms // sites
    n_frames = len(u.trajectory)
    dt       = u.trajectory.dt            # ps

    print(f"  Molecules : {n_mol},  Frames : {n_frames},  dt : {dt} ps")

    # Get average box volume → number density of oxygen atoms
    volumes = []
    for ts in u.trajectory[::max(1, n_frames//200)]:
        volumes.append(ts.volume)          # Å³
    box_vol_nm3 = np.mean(volumes) / 1000.0   # Å³ → nm³
    rho_O = n_mol / box_vol_nm3            # oxygens per nm³

    print(f"  <Box volume> : {box_vol_nm3:.2f} nm³,  ρ_O = {rho_O:.2f} nm⁻³")

    results = {"label": model, "sites": sites, "n_mol": n_mol,
               "n_frames": n_frames, "rho_O_nm3": rho_O}

    # H-atom selection differs between 3-site and 4-site
    h_sel = "name HW1 HW2" if sites == 4 else "name HW1 HW2 or name HW*"

    # ----------------------------------------------------------------------- #
    #  RDF                                                                     #
    # ----------------------------------------------------------------------- #
    print("\n  Computing RDFs ...")
    rdfs = {}
    pairs = [
        ("OO", "name OW", "name OW",   (0.15, 0.80), (1, 1)),
        ("OH", "name OW", h_sel,        (0.05, 0.60), (1, 2) if sites==3 else (1, 2)),
        ("HH", h_sel,      h_sel,        (0.05, 0.60), (2, 2)),
    ]
    for name, s1, s2, rng, excl in pairs:
        try:
            r, g = compute_rdf(u, s1, s2, name, range=rng,
                               exclusion_block=excl)
            rdfs[name] = (r, g)
            print(f"    g_{name}(r): peak = {g.max():.3f} at r = "
                  f"{r[np.argmax(g)]:.3f} nm")
        except Exception as e:
            warnings.warn(f"RDF {name} failed: {e}")

    # ----------------------------------------------------------------------- #
    #  Coordination number                                                     #
    # ----------------------------------------------------------------------- #
    print("\n  Computing coordination numbers ...")
    if "OO" in rdfs:
        r_oo, g_oo = rdfs["OO"]
        cn, r_cut = coordination_number(r_oo, g_oo, rho_O)
        print(f"    CN(O-O) = {cn:.2f}  (r_cut = {r_cut:.3f} nm)")
        results["CN_OO"]    = cn
        results["r_cut_OO"] = r_cut
    if "OH" in rdfs:
        r_oh, g_oh = rdfs["OH"]
        rho_H = 2 * rho_O
        cn_oh, r_cut_oh = coordination_number(r_oh, g_oh, rho_H, r_cut=0.25)
        print(f"    CN(O-H) = {cn_oh:.2f}  (r_cut = {r_cut_oh:.3f} nm)")
        results["CN_OH"] = cn_oh

    # ----------------------------------------------------------------------- #
    #  MSD and diffusion                                                       #
    # ----------------------------------------------------------------------- #
    print("\n  Computing MSD / diffusion ...")
    try:
        times, msd, D = compute_msd(u, sel='name OW')
        print(f"    D = {D:.3e} m²/s  "
              f"(exp water at 300K: ~2.3×10⁻⁹ m²/s)")
        results["D_m2s"]   = D
        results["D_1e9"]   = D * 1e9
    except Exception as e:
        warnings.warn(f"MSD failed: {e}")
        times, msd, D = None, None, None

    # ----------------------------------------------------------------------- #
    #  Tetrahedral order                                                       #
    # ----------------------------------------------------------------------- #
    print("\n  Computing tetrahedral order parameter ...")
    try:
        Q_mean, Q_std = tetrahedral_order(u, n_frames=100)
        print(f"    Q = {Q_mean:.4f} ± {Q_std:.4f}  "
              f"(bulk water ~0.57, ice ~0.95)")
        results["Q_tet_mean"] = Q_mean
        results["Q_tet_std"]  = Q_std
    except Exception as e:
        warnings.warn(f"Tetrahedral order failed: {e}")

    # ----------------------------------------------------------------------- #
    #  Thermodynamic properties from .edr                                      #
    # ----------------------------------------------------------------------- #
    if edr and Path(edr).exists():
        print("\n  Extracting thermodynamic properties ...")
        thermo = extract_energy_properties(edr, out_path, model)
        for name, (mean, std) in thermo.items():
            print(f"    {name:12s}: {mean:.3f} ± {std:.3f}")
            results[f"{name}_mean"] = mean
            results[f"{name}_std"]  = std

    # ----------------------------------------------------------------------- #
    #  Plots                                                                   #
    # ----------------------------------------------------------------------- #
    color = COLORS.get(model.lower(), "#333333")
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle(f"Structural Analysis — {model}  (300 K, 1 bar)",
                 fontsize=14, fontweight="bold")

    # g_OO(r)
    if "OO" in rdfs:
        r, g = rdfs["OO"]
        axes[0,0].plot(r, g, color=color, linewidth=1.8)
        axes[0,0].axvline(results.get("r_cut_OO", 0.35), color="gray",
                          linestyle="--", linewidth=0.8, label="r_cut")
        axes[0,0].set(xlabel="r (nm)", ylabel="g(r)",
                      title="O–O RDF")
        axes[0,0].legend(fontsize=8)
        axes[0,0].grid(True, alpha=0.3)

    # g_OH(r)
    if "OH" in rdfs:
        r, g = rdfs["OH"]
        axes[0,1].plot(r, g, color=color, linewidth=1.8)
        axes[0,1].set(xlabel="r (nm)", ylabel="g(r)", title="O–H RDF")
        axes[0,1].grid(True, alpha=0.3)

    # g_HH(r)
    if "HH" in rdfs:
        r, g = rdfs["HH"]
        axes[0,2].plot(r, g, color=color, linewidth=1.8)
        axes[0,2].set(xlabel="r (nm)", ylabel="g(r)", title="H–H RDF")
        axes[0,2].grid(True, alpha=0.3)

    # Running CN from OO RDF
    if "OO" in rdfs:
        r, g = rdfs["OO"]
        cn_running = cumulative_trapezoid(
            4 * np.pi * rho_O * g * r**2, r, initial=0
        )
        axes[1,0].plot(r, cn_running, color=color, linewidth=1.8)
        axes[1,0].axvline(results.get("r_cut_OO", 0.35), color="gray",
                          linestyle="--", linewidth=0.8)
        axes[1,0].axhline(results.get("CN_OO", 4.5), color="red",
                          linestyle=":", linewidth=0.8,
                          label=f"CN = {results.get('CN_OO', 0):.2f}")
        axes[1,0].set(xlabel="r (nm)", ylabel="CN",
                      title="Running Coordination Number")
        axes[1,0].legend(fontsize=8)
        axes[1,0].grid(True, alpha=0.3)

    # MSD
    if times is not None and msd is not None:
        axes[1,1].plot(times, msd, color=color, linewidth=1.5, alpha=0.8)
        axes[1,1].set(xlabel="t (ps)", ylabel="MSD (nm²)",
                      title=f"MSD  |  D = {D*1e9:.3f} × 10⁻⁹ m²/s")
        axes[1,1].grid(True, alpha=0.3)

    # Thermodynamic summary bar
    thermo_keys = ["Density_mean", "Temperature_mean", "Pressure_mean"]
    thermo_labels = ["Density (kg/m³)", "T (K)", "P (bar)"]
    thermo_vals = [results.get(k, np.nan) for k in thermo_keys]
    valid = [(l, v) for l, v in zip(thermo_labels, thermo_vals)
             if not np.isnan(v)]
    if valid:
        lbls, vals = zip(*valid)
        x = np.arange(len(lbls))
        axes[1,2].bar(x, vals, color=color, alpha=0.8)
        axes[1,2].set_xticks(x)
        axes[1,2].set_xticklabels(lbls, rotation=15, ha="right", fontsize=9)
        axes[1,2].set_title("Thermodynamic Averages")
        axes[1,2].grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    png = out_path / f"analysis_{model}.png"
    plt.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Plot saved → {png}")

    # Save JSON
    # Convert numpy types for JSON serialisation
    def _to_python(obj):
        if isinstance(obj, (np.floating, np.integer)):
            return float(obj)
        return obj

    json_results = {k: _to_python(v) for k, v in results.items()
                    if not isinstance(v, np.ndarray)}
    json_path = out_path / f"results_{model}.json"
    with open(json_path, "w") as fh:
        json.dump(json_results, fh, indent=2)
    print(f"  JSON saved → {json_path}")

    return results, rdfs, (times, msd, D)


# =========================================================================== #
#  Cross-model comparison plot                                                 #
# =========================================================================== #

def plot_comparison(all_results, all_rdfs, all_msds, out_dir="./analysis_results"):
    out_path = Path(out_dir)
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle("Water Model Comparison — 300 K, 1 bar",
                 fontsize=14, fontweight="bold")

    for idx, (res, rdfs, msd_data) in enumerate(zip(all_results, all_rdfs, all_msds)):
        model = res["label"]
        color = COLORS.get(model.lower(), FALLBACK_COLORS[idx % 3])
        lw    = 1.8

        for ax, pair in zip(axes[0], ["OO", "OH", "HH"]):
            if pair in rdfs:
                r, g = rdfs[pair]
                ax.plot(r, g, color=color, linewidth=lw, label=model.upper())

        # Running CN
        if "OO" in rdfs:
            r, g = rdfs["OO"]
            rho  = res["rho_O_nm3"]
            cn_run = cumulative_trapezoid(4*np.pi*rho*g*r**2, r, initial=0)
            axes[1,0].plot(r, cn_run, color=color, linewidth=lw, label=model.upper())

        # MSD
        times, msd, D = msd_data
        if times is not None:
            axes[1,1].plot(times, msd, color=color, linewidth=1.5, alpha=0.9,
                           label=f"{model.upper()}  D={D*1e9:.2f}×10⁻⁹")

    titles = ["O–O RDF", "O–H RDF", "H–H RDF",
              "Coordination Number", "MSD", "Thermodynamics"]
    xlabels = ["r (nm)"]*3 + ["r (nm)", "t (ps)", "Model"]
    ylabels = ["g(r)"]*3 + ["CN", "MSD (nm²)", ""]

    for ax, title, xl, yl in zip(axes.flat, titles, xlabels, ylabels):
        ax.set(title=title, xlabel=xl, ylabel=yl)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    # Thermodynamic bar chart
    ax_thermo = axes[1, 2]
    keys    = ["Density_mean", "Temperature_mean"]
    ylabels_t = ["Density (kg/m³)", "T (K)"]
    x = np.arange(len(all_results))
    width = 0.35
    for ki, (key, ylabel_t) in enumerate(zip(keys, ylabels_t)):
        vals = [r.get(key, np.nan) for r in all_results]
        labels = [r["label"].upper() for r in all_results]
        if not all(np.isnan(vals)):
            ax2 = ax_thermo if ki == 0 else ax_thermo.twinx()
            clrs = [COLORS.get(r["label"].lower(), FALLBACK_COLORS[i])
                    for i, r in enumerate(all_results)]
            ax2.bar(x + ki*width, vals, width, color=clrs, alpha=0.7,
                    label=ylabel_t)
            ax2.set_ylabel(ylabel_t, fontsize=8)

    ax_thermo.set_xticks(x + width/2)
    ax_thermo.set_xticklabels([r["label"].upper() for r in all_results])
    ax_thermo.set_title("Thermodynamic Comparison")

    plt.tight_layout()
    png = out_path / "comparison_all_models.png"
    plt.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nComparison plot saved → {png}")


# =========================================================================== #
#  CLI                                                                         #
# =========================================================================== #

def parse_args():
    p = argparse.ArgumentParser(
        description="Full structural analysis of GROMACS water trajectories"
    )
    p.add_argument("--model",   help="Model label (e.g. spcfw)")
    p.add_argument("--traj",    help="Path to .xtc trajectory")
    p.add_argument("--tpr",     help="Path to .tpr file")
    p.add_argument("--edr",     help="Path to .edr energy file")
    p.add_argument("--sites",   type=int, default=3,
                   help="Atoms per molecule: 3 (SPC/Fw, FBA/eps) or 4 (TIP4P/2005f)")
    p.add_argument("--out_dir", default="./analysis_results")
    # Comparison mode
    p.add_argument("--compare", action="store_true",
                   help="Compare multiple models")
    p.add_argument("--models",  nargs="+",
                   help="Model labels for comparison")
    p.add_argument("--dirs",    nargs="+",
                   help="npt/ directory for each model (contains .xtc, .tpr, .edr)")
    p.add_argument("--sites_list", nargs="+", type=int,
                   help="Sites per molecule for each model")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if not HAS_MDA:
        print("ERROR: MDAnalysis is required.")
        print("Install with:  pip install MDAnalysis")
        import sys; sys.exit(1)

    if args.compare:
        # Compare mode
        models  = args.models  or ["spcfw", "fbaeps", "tip4p2005f"]
        dirs    = args.dirs    or [f"retis_{m}/npt" for m in models]
        sites_l = args.sites_list or [3, 3, 4]

        all_results, all_rdfs, all_msds = [], [], []
        for model, d, s in zip(models, dirs, sites_l):
            dp   = Path(d)
            traj = next(dp.glob("*.xtc"), None)
            tpr  = next(dp.glob("*.tpr"), None)
            edr  = next(dp.glob("*.edr"), None)
            if traj is None or tpr is None:
                print(f"WARNING: missing traj/tpr in {d}, skipping {model}")
                continue
            res, rdfs, msd_data = analyse_model(
                model, traj, tpr, edr, s, args.out_dir
            )
            all_results.append(res)
            all_rdfs.append(rdfs)
            all_msds.append(msd_data)

        if len(all_results) > 1:
            plot_comparison(all_results, all_rdfs, all_msds, args.out_dir)

    else:
        # Single model mode
        if not args.model or not args.traj or not args.tpr:
            print("Provide --model, --traj, --tpr  (and optionally --edr)")
            print("Or use --compare mode.")
            import sys; sys.exit(1)

        analyse_model(args.model, args.traj, args.tpr, args.edr,
                      args.sites, args.out_dir)
