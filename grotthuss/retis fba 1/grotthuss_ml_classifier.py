"""
Grotthuss/Autoionization Initiation Classifier
Based on: Moqadam et al., PNAS 2018

This script:
1. Reads accepted PyRETIS trajectory frames
2. Computes the 4 key CVs from the paper:
     w4   - 4-water hydrogen bond wire length (Å)
     q    - tetrahedral order parameter of Oλ
     na   - number of H-bonds accepted by Oλ
     qcos - angle order parameter of the wire
3. Applies the decision tree rules from Fig. 5 of the paper
4. Classifies each frame as "reactive-likely" or "non-reactive"
5. Plots distributions matching Fig. 4 of the paper

Run from PyRETIS simulation folder:
    python3 grotthuss_ml_classifier.py

Requires: MDAnalysis, numpy, matplotlib, scipy, scikit-learn
    pip install scikit-learn --break-system-packages
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.spatial.distance import cdist
import warnings, os, glob
warnings.filterwarnings('ignore')

import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array, calc_angles

# ── Configuration ─────────────────────────────────────────────────────────────
GRO          = 'gromacs_input/conf.gro'   # structure file
TRAJ_PATTERN = '*/traj/traj-acc/*.trr'    # all accepted path trajectories
LAMBDA_C     = 0.109                       # nm — initiation threshold from paper (1.15 Å)
N_ENSEMBLES  = 7                           # number of PyRETIS path ensembles

# Paper's decision tree thresholds (Fig. 5A, converted to nm where needed):
# w4 < 7.6 Å = 0.76 nm → likely reactive (compressed wire)
# na = 4 → strongly reactive (hypercoordinated)
# q < 0.5 → distorted tetrahedral → initiates transfer
W4_THRESHOLD   = 0.76    # nm (7.6 Å from paper)
W4_COMPRESSED  = 0.73    # nm (7.3 Å from paper — "strongly compressed")
NA_HYPER       = 4       # accepted H-bonds = hypercoordination
Q_DISTORTED    = 0.5     # below this = distorted tetrahedral

# H-bond geometric criteria (Luzar & Chandler 1996, used in paper)
HB_OO_MAX   = 0.35    # nm — O-O distance cutoff
HB_OHO_MIN  = 150.0   # degrees — O-H...O angle cutoff

print("=" * 65)
print("Grotthuss Initiation Classifier — Moqadam et al. 2018")
print("=" * 65)

# ── Load structure ─────────────────────────────────────────────────────────────
if not os.path.exists(GRO):
    raise FileNotFoundError(f"Structure file not found: {GRO}")

u_ref = mda.Universe(GRO)
n_atoms = len(u_ref.atoms)
n_mol   = n_atoms // 3
print(f"System: {n_mol} water molecules, {n_atoms} atoms")

# ── Find trajectory files ──────────────────────────────────────────────────────
traj_files = sorted(glob.glob(TRAJ_PATTERN))
if not traj_files:
    # fallback: look for any trr in ensemble folders
    traj_files = sorted(glob.glob('00*/traj/traj-acc/*.trr'))
print(f"Found {len(traj_files)} accepted trajectory segments")

# ── CV computation functions ───────────────────────────────────────────────────

def find_lambda_oxygen(ow_pos, hw_pos, box):
    """
    Find Oλ: the oxygen with the largest covalent O-H distance.
    Returns (o_lambda_idx, max_oh_dist_nm)
    Paper: "oxygen atom for which the covalent O–H distance is largest"
    """
    max_oh = 0.0
    o_lambda = 0
    n_o = len(ow_pos)
    for i in range(n_o):
        h1 = hw_pos[2*i]
        h2 = hw_pos[2*i + 1]
        d1 = h1 - ow_pos[i]
        d2 = h2 - ow_pos[i]
        # minimum image
        for d in range(3):
            L = box[d]
            d1[d] -= L * np.round(d1[d] / L)
            d2[d] -= L * np.round(d2[d] / L)
        r1 = np.linalg.norm(d1)
        r2 = np.linalg.norm(d2)
        r_max = max(r1, r2)
        if r_max > max_oh:
            max_oh = r_max
            o_lambda = i
    return o_lambda, max_oh


def compute_hbonds(ow_pos, hw_pos, box):
    """
    Compute hydrogen bond network using Luzar-Chandler criteria.
    Returns adjacency list: hbonds[i] = list of j where i donates to j
    and na[i] = number of H-bonds accepted by molecule i
    """
    n_mol = len(ow_pos)
    # O-O distance matrix
    oo_dist = distance_array(ow_pos, ow_pos, box=box)
    np.fill_diagonal(oo_dist, np.inf)

    # acceptors[i] = list of molecules that donate to i
    acceptors = [[] for _ in range(n_mol)]
    donors    = [[] for _ in range(n_mol)]

    for i in range(n_mol):
        # potential acceptors: O-O < 3.5 Å
        candidates = np.where(oo_dist[i] < HB_OO_MAX)[0]
        for j in candidates:
            # check both H atoms of molecule i
            for h_local in range(2):
                h_pos = hw_pos[2*i + h_local]
                o_i   = ow_pos[i]
                o_j   = ow_pos[j]
                # O-H...O angle
                vec_oh = h_pos - o_i
                vec_ha = o_j - h_pos
                for d in range(3):
                    L = box[d]
                    vec_oh[d] -= L * np.round(vec_oh[d] / L)
                    vec_ha[d] -= L * np.round(vec_ha[d] / L)
                cos_a = np.dot(vec_oh, vec_ha) / (
                    np.linalg.norm(vec_oh) * np.linalg.norm(vec_ha) + 1e-10)
                angle = np.degrees(np.arccos(np.clip(cos_a, -1, 1)))
                if angle > HB_OHO_MIN:
                    donors[i].append(j)
                    acceptors[j].append(i)

    na = np.array([len(acceptors[i]) for i in range(n_mol)])
    nd = np.array([len(donors[i]) for i in range(n_mol)])
    return donors, acceptors, na, nd


def find_wire_length(o_lambda, ow_pos, hb_donors, box, wire_size=4):
    """
    Find shortest H-bond wire of length wire_size containing Oλ.
    Returns wire length w_i (sum of O-O distances of consecutive members) in nm.
    Paper: w4 = length of 4-water wire containing Oλ
    """
    n_mol = len(ow_pos)
    # BFS to find all paths of length wire_size-1 from o_lambda
    best_length = np.inf
    # Simple BFS up to depth wire_size-1
    queue = [(o_lambda, [o_lambda])]
    while queue:
        node, path = queue.pop(0)
        if len(path) == wire_size:
            # compute wire length
            length = 0.0
            for k in range(len(path)-1):
                d = ow_pos[path[k+1]] - ow_pos[path[k]]
                for dim in range(3):
                    L = box[dim]
                    d[dim] -= L * np.round(d[dim] / L)
                length += np.linalg.norm(d)
            best_length = min(best_length, length)
            continue
        for neighbor in hb_donors[node]:
            if neighbor not in path:
                queue.append((neighbor, path + [neighbor]))
    return best_length if best_length < np.inf else -1.0


def compute_tetrahedral_q(o_lambda, ow_pos, box):
    """
    Tetrahedral order parameter q of Oλ.
    q = 1 - (3/8) Σ_{j<k} (cos ψ_jk + 1/3)²
    where ψ_jk is the angle formed by Oλ and its 4 nearest oxygen neighbors.
    Paper Eq. 3; q=1 for perfect tetrahedron.
    """
    D = distance_array(ow_pos[o_lambda:o_lambda+1], ow_pos, box=box)[0]
    D[o_lambda] = np.inf
    nn_idx = np.argsort(D)[:4]
    v = ow_pos[nn_idx] - ow_pos[o_lambda]
    for dim in range(3):
        L = box[dim]
        v[:, dim] -= L * np.round(v[:, dim] / L)
    norms = np.linalg.norm(v, axis=1, keepdims=True).clip(1e-10)
    v /= norms
    s = sum((np.dot(v[j], v[k]) + 1/3)**2
            for j in range(3) for k in range(j+1, 4))
    return 1.0 - (3/8) * s


def compute_qcos(o_lambda, wire_oxygens, ow_pos, box):
    """
    Angular order parameter qcos = min(cos α, cos β)
    where α and β are the two internal angles in the wire.
    Paper: "smallest of the cosine of the two internal angles in the wire"
    """
    if len(wire_oxygens) < 3:
        return -1.0
    # Use first 3 members of wire for internal angles
    w = [ow_pos[i].copy() for i in wire_oxygens[:3]]
    # Angle at middle atom
    v1 = w[0] - w[1]
    v2 = w[2] - w[1]
    for dim in range(3):
        L = box[dim]
        v1[dim] -= L * np.round(v1[dim] / L)
        v2[dim] -= L * np.round(v2[dim] / L)
    cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
    # If wire has 4 members, also compute second internal angle
    if len(wire_oxygens) >= 4:
        w2 = [ow_pos[i].copy() for i in wire_oxygens[1:4]]
        v3 = w2[0] - w2[1]
        v4 = w2[2] - w2[1]
        for dim in range(3):
            L = box[dim]
            v3[dim] -= L * np.round(v3[dim] / L)
            v4[dim] -= L * np.round(v4[dim] / L)
        cos_b = np.dot(v3, v4) / (np.linalg.norm(v3) * np.linalg.norm(v4) + 1e-10)
        return min(cos_a, cos_b)
    return cos_a


def apply_decision_tree(w4, q, na, qcos, lambda_val):
    """
    Decision tree from Fig. 5A of Moqadam et al. 2018.

    Key rules (simplified from paper's CART tree):
    Rule 1: w4 >= 7.6 Å → NON-REACTIVE (wire too long)
    Rule 2: w4 < 7.6 Å AND λ2 < 1.1 Å → check further
    Rule 3: w4 < 7.3 Å AND na >= 2.5 → REACTIVE-LIKELY
    Rule 4: w4 < 7.3 Å AND na >= 4 AND q4 <= 0.23 → STRONGLY REACTIVE
    Rule 5: w4 >= 7.6 Å AND λ2 >= 1.1 Å → REACTIVE (rule 5 in paper)

    Returns: (label, probability_estimate)
        label = 'reactive' | 'nonreactive' | 'uncertain'
        prob  = estimated probability of being reactive
    """
    w4_A = w4 * 10  # convert nm to Å for comparison with paper thresholds

    # Rule 1: wire too long → non-reactive
    if w4_A >= 7.6 and lambda_val * 10 < 1.1:
        return 'nonreactive', 0.13   # 13/87 from paper node 3

    # Rule 5: large wire but stretched bond → reactive
    if w4_A >= 7.6 and lambda_val * 10 >= 1.1:
        return 'reactive', 0.75      # node 2: 75/25

    # Compressed wire region (w4 < 7.6 Å)
    if w4_A < 7.3 and na >= NA_HYPER:
        return 'reactive', 0.15      # node: 15% chance — paper's "best" condition

    if w4_A < 7.5 and na >= 3:
        return 'reactive', 0.04      # moderate

    if w4_A < 7.6 and q < Q_DISTORTED:
        return 'uncertain', 0.02

    return 'nonreactive', 0.005


# ── Main analysis loop ─────────────────────────────────────────────────────────
results = []  # list of dicts, one per frame analyzed

print(f"\nAnalyzing frames near λ_c = {LAMBDA_C*10:.2f} Å...")
n_frames_analyzed = 0

for traj_file in traj_files[:50]:   # limit to first 50 files for speed
    try:
        u = mda.Universe(GRO, traj_file)
        OW = u.select_atoms('name OW')
        HW = u.select_atoms('name HW1 or name HW2')

        for ts in u.trajectory:
            ow_pos = OW.positions / 10.0   # Å → nm
            hw_pos = HW.positions / 10.0
            box    = ts.dimensions[:3] / 10.0  # Å → nm

            # Find Oλ and λ value
            o_lambda, lambda_val = find_lambda_oxygen(ow_pos, hw_pos, box)

            # Only analyze frames near the initiation threshold
            if abs(lambda_val - LAMBDA_C) > 0.010:
                continue

            # Compute H-bond network
            donors, acceptors, na_arr, nd_arr = compute_hbonds(
                ow_pos, hw_pos, box)

            na_lambda = na_arr[o_lambda]

            # Find wire containing Oλ
            w4 = find_wire_length(o_lambda, ow_pos, donors, box, wire_size=4)

            # Tetrahedral order of Oλ
            q_lambda = compute_tetrahedral_q(o_lambda, ow_pos, box)

            # Wire angle parameter (simplified)
            qcos_lambda = -1.0   # computed below if wire found
            if w4 > 0:
                # Get actual wire members for qcos
                wire = [o_lambda]
                visited = {o_lambda}
                for _ in range(3):
                    if not donors[wire[-1]]:
                        break
                    for nxt in donors[wire[-1]]:
                        if nxt not in visited:
                            wire.append(nxt)
                            visited.add(nxt)
                            break
                qcos_lambda = compute_qcos(o_lambda, wire, ow_pos, box)

            # Apply decision tree
            label, prob = apply_decision_tree(
                w4, q_lambda, na_lambda, qcos_lambda, lambda_val)

            results.append({
                'frame':    ts.frame,
                'traj':     os.path.basename(traj_file),
                'lambda':   lambda_val * 10,    # Å for display
                'o_lambda': o_lambda,
                'w4':       w4 * 10,            # nm → Å
                'q':        q_lambda,
                'na':       na_lambda,
                'qcos':     qcos_lambda,
                'label':    label,
                'prob':     prob,
            })
            n_frames_analyzed += 1

    except Exception as e:
        print(f"  Warning: could not process {traj_file}: {e}")
        continue

print(f"\nAnalyzed {n_frames_analyzed} frames near λ_c")

if not results:
    print("\nNo frames found near λ_c threshold.")
    print("This is expected if O-H bonds don't stretch that far with SPC/Fw.")
    print("Try lowering LAMBDA_C to e.g. 0.110 nm (1.10 Å)")
    print("Or run the standalone analysis below on the equilibration trajectory.")
    # ── Fallback: analyze production trajectory directly ──────────────────────
    print("\nRunning fallback analysis on production trajectory...")
    prod_gro = 'equil/npt_300K.gro' if os.path.exists('equil/npt_300K.gro') \
               else '../equil/npt_300K.gro'
    prod_xtc = 'equil/npt_300K.xtc' if os.path.exists('equil/npt_300K.xtc') \
               else '../equil/npt_300K.xtc'

    if os.path.exists(prod_gro) and os.path.exists(prod_xtc):
        u = mda.Universe(prod_gro, prod_xtc)
        OW = u.select_atoms('name OW')
        HW = u.select_atoms('name HW1 or name HW2')

        lambda_vals = []
        for ts in u.trajectory:
            ow_pos = OW.positions / 10.0
            hw_pos = HW.positions / 10.0
            box    = ts.dimensions[:3] / 10.0
            _, lv  = find_lambda_oxygen(ow_pos, hw_pos, box)
            lambda_vals.append(lv * 10)  # Å

        lambda_vals = np.array(lambda_vals)
        print(f"\nλ statistics from production trajectory:")
        print(f"  Mean:     {np.mean(lambda_vals):.4f} Å")
        print(f"  Std:      {np.std(lambda_vals):.4f} Å")
        print(f"  Max:      {np.max(lambda_vals):.4f} Å")
        print(f"  Min:      {np.min(lambda_vals):.4f} Å")
        print(f"  >1.10 Å:  {np.sum(lambda_vals > 1.10)} frames ({100*np.mean(lambda_vals > 1.10):.2f}%)")
        print(f"  >1.13 Å:  {np.sum(lambda_vals > 1.13)} frames ({100*np.mean(lambda_vals > 1.13):.2f}%)")
        print(f"  >1.15 Å:  {np.sum(lambda_vals > 1.15)} frames")
        print(f"\nSuggested interfaces for SPC/Fw based on distribution:")
        pcts = [50, 75, 90, 95, 99]
        for p in pcts:
            print(f"  {p}th percentile: {np.percentile(lambda_vals, p):.4f} Å")

        # Plot λ distribution
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(lambda_vals, bins=100, density=True,
                color='steelblue', edgecolor='white', lw=0.3)
        ax.axvline(1.10, color='red', ls='--', lw=1.5, label='1.10 Å')
        ax.axvline(1.13, color='orange', ls='--', lw=1.5, label='1.13 Å')
        ax.axvline(1.15, color='green', ls='--', lw=1.5,
                   label='1.15 Å (paper threshold)')
        ax.set_xlabel('λ (Å) — largest O-H distance', fontsize=12)
        ax.set_ylabel('P(λ)', fontsize=12)
        ax.set_title('Distribution of λ in SPC/Fw water\n'
                     'Dashed lines = suggested PyRETIS interfaces', fontsize=11)
        ax.legend()
        plt.tight_layout()
        plt.savefig('lambda_distribution.png', dpi=150)
        print("\nSaved: lambda_distribution.png")
        print("→ Use this to set correct interfaces in retis_grotthuss.rst")
    else:
        print(f"Production trajectory not found at {prod_xtc}")

else:
    # ── Results summary ────────────────────────────────────────────────────────
    reactive     = [r for r in results if r['label'] == 'reactive']
    nonreactive  = [r for r in results if r['label'] == 'nonreactive']
    uncertain    = [r for r in results if r['label'] == 'uncertain']

    print(f"\nClassification results:")
    print(f"  Reactive-likely:  {len(reactive):4d} ({100*len(reactive)/len(results):.1f}%)")
    print(f"  Uncertain:        {len(uncertain):4d} ({100*len(uncertain)/len(results):.1f}%)")
    print(f"  Non-reactive:     {len(nonreactive):4d} ({100*len(nonreactive)/len(results):.1f}%)")

    # ── Plot: Fig 4 equivalent ─────────────────────────────────────────────────
    fig = plt.figure(figsize=(16, 10))
    fig.suptitle('Water Autoionization Initiation — Moqadam et al. 2018\n'
                 f'SPC/Fw, {len(results)} frames analyzed near λ_c={LAMBDA_C*10:.2f} Å',
                 fontsize=13)
    gs = gridspec.GridSpec(2, 3, hspace=0.4, wspace=0.35)

    def split_cv(key):
        r_vals = [r[key] for r in reactive   if r[key] > -1]
        n_vals = [r[key] for r in nonreactive if r[key] > -1]
        return np.array(r_vals), np.array(n_vals)

    # w4 distribution (main discriminator, Fig 4A)
    ax1 = fig.add_subplot(gs[0, 0])
    rw, nw = split_cv('w4')
    if len(rw): ax1.hist(rw, bins=20, density=True, alpha=0.6,
                         color='blue', label='Reactive', edgecolor='white')
    if len(nw): ax1.hist(nw, bins=20, density=True, alpha=0.6,
                         color='gray', label='Non-reactive', edgecolor='white')
    ax1.axvline(7.6, color='red', ls='--', lw=1.5, label='7.6 Å threshold')
    ax1.axvline(7.3, color='orange', ls=':', lw=1.5, label='7.3 Å (strongly comp.)')
    ax1.set_xlabel('w₄ (Å)'); ax1.set_ylabel('P(w₄)')
    ax1.set_title('Wire length w₄\n(most important — Fig 3A)'); ax1.legend(fontsize=7)

    # na distribution
    ax2 = fig.add_subplot(gs[0, 1])
    rn, nn = split_cv('na')
    bins_na = np.arange(-0.5, 7.5, 1)
    if len(rn): ax2.hist(rn, bins=bins_na, density=True, alpha=0.6,
                         color='blue', label='Reactive', edgecolor='white')
    if len(nn): ax2.hist(nn, bins=bins_na, density=True, alpha=0.6,
                         color='gray', label='Non-reactive', edgecolor='white')
    ax2.axvline(3.5, color='red', ls='--', lw=1.5, label='hypercoord. (na=4)')
    ax2.set_xlabel('na (H-bonds accepted by Oλ)'); ax2.set_ylabel('P(na)')
    ax2.set_title('H-bonds accepted\n(Fig 4A)'); ax2.legend(fontsize=7)

    # q distribution
    ax3 = fig.add_subplot(gs[0, 2])
    rq, nq = split_cv('q')
    if len(rq): ax3.hist(rq, bins=20, density=True, alpha=0.6,
                         color='blue', label='Reactive', edgecolor='white')
    if len(nq): ax3.hist(nq, bins=20, density=True, alpha=0.6,
                         color='gray', label='Non-reactive', edgecolor='white')
    ax3.axvline(0.5, color='red', ls='--', lw=1.5, label='q=0.5')
    ax3.set_xlabel('q (tetrahedral OP of Oλ)'); ax3.set_ylabel('P(q)')
    ax3.set_title('Tetrahedral order\n(Fig 3A)'); ax3.legend(fontsize=7)

    # 2D: w4 vs na (Fig 4A reproduction)
    ax4 = fig.add_subplot(gs[1, 0:2])
    rw4   = [r['w4'] for r in reactive]
    rna   = [r['na'] for r in reactive]
    nw4   = [r['w4'] for r in nonreactive]
    nna   = [r['na'] for r in nonreactive]
    if nw4: ax4.scatter(nw4, nna, c='gray',  alpha=0.3, s=15, label='Non-reactive')
    if rw4: ax4.scatter(rw4, rna, c='blue',  alpha=0.6, s=25, label='Reactive')
    ax4.axvline(7.6, color='red', ls='--', lw=1.5)
    ax4.axvline(7.3, color='orange', ls=':', lw=1.5)
    ax4.axhline(3.5, color='green', ls=':', lw=1.5)
    ax4.set_xlabel('w₄ (Å)'); ax4.set_ylabel('na')
    ax4.set_title('w₄ vs na — Fig 4A reproduction\nBlue=reactive, Gray=non-reactive')
    ax4.legend(fontsize=8)

    # Summary
    ax5 = fig.add_subplot(gs[1, 2])
    ax5.axis('off')
    if reactive:
        best = max(reactive, key=lambda r: r['prob'])
        summary = (
            f"SUMMARY\n{'─'*35}\n"
            f"Frames analyzed: {len(results)}\n"
            f"Reactive-likely: {len(reactive)} ({100*len(reactive)/len(results):.0f}%)\n"
            f"Non-reactive:    {len(nonreactive)} ({100*len(nonreactive)/len(results):.0f}%)\n\n"
            f"BEST REACTIVE FRAME:\n"
            f"  Traj:   {best['traj']}\n"
            f"  λ:      {best['lambda']:.4f} Å\n"
            f"  w₄:     {best['w4']:.3f} Å\n"
            f"  na:     {best['na']}\n"
            f"  q:      {best['q']:.4f}\n"
            f"  P(rxn): ~{best['prob']:.3f}\n\n"
            f"Paper: P(rxn|w4<7.3,na=4) = 0.15\n"
            f"Paper: P(rxn|random) = 10⁻⁷"
        )
    else:
        summary = "No reactive frames found.\nLower LAMBDA_C or run longer."
    ax5.text(0.05, 0.95, summary, transform=ax5.transAxes, fontsize=9,
             va='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.savefig('grotthuss_analysis.png', dpi=150, bbox_inches='tight')
    print("\nSaved: grotthuss_analysis.png")

    # Save results to CSV
    import csv
    with open('grotthuss_results.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print("Saved: grotthuss_results.csv")

print("\nDone.")
print("\nNext steps:")
print("  1. Check lambda_distribution.png to set correct interfaces")
print("     for retis_grotthuss.rst")
print("  2. Run PyRETIS: pyretisrun -i retis_grotthuss.rst -p")
print("  3. Run this script again on the PyRETIS output")
print("  4. Frames labeled 'reactive' are Grotthuss-initiation candidates")
