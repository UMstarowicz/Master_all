#!/usr/bin/env python3
"""
O-H distance analysis for SPC/Fw water RETIS simulation.
Extracts O-H distances from trajectory, computes P(r) and free energy G(r).

Usage:
    python analyze_oh.py --trr path/to/traj.trr --gro path/to/conf.gro
    python analyze_oh.py --trr load/prd_reduced_2.trr --gro gromacs_input/conf.gro

Output:
    oh_analysis.png  — histogram, log(P), free energy plots
    oh_distances.txt — raw max O-H distances per frame
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys

# ── Parse arguments ────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description='O-H distance analysis')
parser.add_argument('--trr', default='load/prd_reduced_2.trr',
                    help='GROMACS trajectory file (.trr or .xtc)')
parser.add_argument('--gro', default='gromacs_input/conf.gro',
                    help='Reference structure file (.gro)')
parser.add_argument('--temp', type=float, default=270.0,
                    help='Temperature in K (default: 270)')
parser.add_argument('--out', default='oh_analysis.png',
                    help='Output plot filename')
parser.add_argument('--cutoff', type=float, default=0.14,
                    help='Max covalent O-H cutoff in nm (default: 0.14)')
args = parser.parse_args()

# ── Check files ────────────────────────────────────────────────────────────────
for f in [args.trr, args.gro]:
    if not os.path.exists(f):
        print(f"ERROR: File not found: {f}")
        sys.exit(1)

print(f"Loading trajectory: {args.trr}")
print(f"Reference structure: {args.gro}")

# ── Load with MDTraj ───────────────────────────────────────────────────────────
try:
    import mdtraj as md
except ImportError:
    print("ERROR: mdtraj not installed. Run: pip install mdtraj")
    sys.exit(1)

top  = md.load_topology(args.gro)
traj = md.load(args.trr, top=args.gro)
print(f"Loaded {traj.n_frames} frames, {traj.n_atoms} atoms")

# ── Extract O-H distances ──────────────────────────────────────────────────────
# SPC/Fw atom order per molecule: OW(3i), HW1(3i+1), HW2(3i+2)
# Units: nm (MDTraj native)

n_mol   = traj.n_atoms // 3
pos     = traj.xyz  # shape: (n_frames, n_atoms, 3)

print(f"Computing O-H distances for {n_mol} molecules × {traj.n_frames} frames...")

# Build O-H1 and O-H2 pairs
oh1_pairs = np.array([[3*i, 3*i+1] for i in range(n_mol)])
oh2_pairs = np.array([[3*i, 3*i+2] for i in range(n_mol)])

d_oh1 = md.compute_distances(traj, oh1_pairs)  # (n_frames, n_mol)
d_oh2 = md.compute_distances(traj, oh2_pairs)  # (n_frames, n_mol)

# All O-H distances (every bond, every frame)
all_oh = np.concatenate([d_oh1.ravel(), d_oh2.ravel()])

# Max O-H per frame (the order parameter)
max_oh_per_frame = np.maximum(d_oh1.max(axis=1), d_oh2.max(axis=1))

# Filter to covalent range only
covalent_mask = all_oh < args.cutoff
oh_covalent   = all_oh[covalent_mask]

print(f"Total O-H bonds analysed: {len(all_oh):,}")
print(f"Covalent O-H bonds (< {args.cutoff} nm): {len(oh_covalent):,}")
print(f"\nO-H distance statistics (nm):")
print(f"  Mean:   {oh_covalent.mean()*10:.4f} Å")
print(f"  Std:    {oh_covalent.std()*10:.4f} Å")
print(f"  Max:    {oh_covalent.max()*10:.4f} Å  ({oh_covalent.max():.5f} nm)")
print(f"\nMax O-H per frame:")
print(f"  Mean:   {max_oh_per_frame.mean()*10:.4f} Å")
print(f"  Max:    {max_oh_per_frame.max()*10:.4f} Å  ({max_oh_per_frame.max():.5f} nm)")

# Save raw data
np.savetxt('oh_distances.txt',
           np.column_stack([np.arange(traj.n_frames), max_oh_per_frame]),
           header='frame  max_OH_nm', fmt=['%d', '%.8f'])
print(f"\nSaved per-frame max O-H to: oh_distances.txt")

# ── Histogram and free energy ──────────────────────────────────────────────────
kB   = 0.008314  # kJ/mol/K
T    = args.temp
kBT  = kB * T

# Use all covalent O-H distances for the distribution
bins     = np.linspace(0.095, args.cutoff, 200)
counts, edges = np.histogram(oh_covalent, bins=bins)
centers  = 0.5 * (edges[:-1] + edges[1:])

# Probability density (normalized)
dr       = edges[1] - edges[0]
prob     = counts / (counts.sum() * dr)

# Free energy G(r) = -kBT * ln(P(r)), set minimum to zero
with np.errstate(divide='ignore', invalid='ignore'):
    log_prob = np.where(prob > 0, np.log(prob), np.nan)
    G        = -kBT * log_prob
    G       -= np.nanmin(G)

# Theoretical harmonic free energy for SPC/Fw
r0   = 0.1012   # nm
kb   = 886062.0 # kJ/mol/nm^2
r_th = np.linspace(0.095, 0.135, 300)
G_th = 0.5 * kb * (r_th - r0)**2  # harmonic PE
G_th -= G_th.min()

# ── Plot ───────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle(f'SPC/Fw O–H Bond Analysis  |  T = {T:.0f} K  |  {traj.n_frames} frames',
             fontsize=14, fontweight='bold', y=0.98)

colors = {'hist': '#2196F3', 'logp': '#E91E63', 'free': '#4CAF50',
          'harm': '#FF9800', 'op': '#9C27B0'}

# ── Panel 1: Histogram ─────────────────────────────────────────────────────────
ax = axes[0, 0]
ax.bar(centers * 10, counts, width=dr * 10, color=colors['hist'],
       alpha=0.75, edgecolor='none', label='O–H bonds')
ax.axvline(r0 * 10, color='k', ls='--', lw=1.5, label=f'r₀ = {r0*10:.2f} Å')
ax.axvline(args.cutoff * 10, color='red', ls=':', lw=1.5,
           label=f'cutoff = {args.cutoff*10:.1f} Å')
ax.set_xlabel('O–H distance (Å)', fontsize=11)
ax.set_ylabel('Count', fontsize=11)
ax.set_title('O–H Distance Histogram', fontsize=12)
ax.legend(fontsize=9)
ax.set_xlim(0.95, args.cutoff * 10 + 0.1)

# ── Panel 2: log P(r) ─────────────────────────────────────────────────────────
ax = axes[0, 1]
valid = np.isfinite(log_prob)
ax.plot(centers[valid] * 10, log_prob[valid], color=colors['logp'],
        lw=2, label='log P(r)')
ax.axvline(r0 * 10, color='k', ls='--', lw=1.5, label=f'r₀ = {r0*10:.2f} Å')
# Mark the max sampled distance
max_sampled = oh_covalent.max()
ax.axvline(max_sampled * 10, color='red', ls='-', lw=2,
           label=f'Max sampled = {max_sampled*10:.3f} Å')
ax.set_xlabel('O–H distance (Å)', fontsize=11)
ax.set_ylabel('log P(r)', fontsize=11)
ax.set_title('Log Probability', fontsize=12)
ax.legend(fontsize=9)
ax.set_xlim(0.95, args.cutoff * 10 + 0.1)

# ── Panel 3: Free energy G(r) ─────────────────────────────────────────────────
ax = axes[1, 0]
valid = np.isfinite(G)
ax.plot(centers[valid] * 10, G[valid], color=colors['free'],
        lw=2.5, label='G(r) = −kBT ln P(r)')
ax.plot(r_th * 10, G_th, color=colors['harm'], lw=1.5,
        ls='--', alpha=0.8, label='Harmonic PE (SPC/Fw)')
ax.axhline(kBT, color='grey', ls=':', lw=1, alpha=0.7, label=f'kT = {kBT:.2f} kJ/mol')
ax.axhline(2*kBT, color='grey', ls=':', lw=1, alpha=0.5)
ax.axhline(3*kBT, color='grey', ls=':', lw=1, alpha=0.3)
ax.set_xlabel('O–H distance (Å)', fontsize=11)
ax.set_ylabel('G(r)  (kJ/mol)', fontsize=11)
ax.set_title('Free Energy Profile', fontsize=12)
ax.legend(fontsize=9)
ax.set_xlim(0.95, args.cutoff * 10 + 0.1)
ax.set_ylim(-1, 20 * kBT)

# Add kT labels
for n, label in [(1, '1 kT'), (2, '2 kT'), (3, '3 kT')]:
    ax.text(args.cutoff * 10 + 0.05, n * kBT, label,
            va='center', fontsize=8, color='grey')

# ── Panel 4: Max O-H order parameter over time ────────────────────────────────
ax = axes[1, 1]
frames = np.arange(traj.n_frames)
ax.plot(frames, max_oh_per_frame * 10, color=colors['op'],
        lw=0.8, alpha=0.8, label='max O–H (order param)')
ax.axhline(r0 * 10, color='k', ls='--', lw=1.5, label=f'r₀ = {r0*10:.2f} Å')
# Draw interface lines
interfaces_nm = [0.1085, 0.1090, 0.1095, 0.1100, 0.1105,
                 0.1110, 0.1115, 0.1120]
for iface in interfaces_nm:
    ax.axhline(iface * 10, color='#2196F3', ls=':', lw=0.7, alpha=0.5)
ax.axhline(interfaces_nm[0] * 10, color='#2196F3', ls=':', lw=0.7,
           alpha=0.5, label='RETIS interfaces')
ax.set_xlabel('Frame', fontsize=11)
ax.set_ylabel('Max O–H distance (Å)', fontsize=11)
ax.set_title('Order Parameter over Trajectory', fontsize=12)
ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig(args.out, dpi=150, bbox_inches='tight')
print(f"Saved plot to: {args.out}")

# ── Summary ────────────────────────────────────────────────────────────────────
print(f"\n{'='*55}")
print(f"SUMMARY")
print(f"{'='*55}")
print(f"Temperature:          {T:.0f} K")
print(f"kBT:                  {kBT:.4f} kJ/mol")
print(f"Equilibrium r₀:       {r0*10:.3f} Å")
print(f"Max sampled O-H:      {oh_covalent.max()*10:.4f} Å")
print(f"Max OP per frame:     {max_oh_per_frame.max()*10:.4f} Å")

# Energy at max sampled
r_max = oh_covalent.max()
E_max = 0.5 * kb * (r_max - r0)**2
print(f"Energy at max O-H:    {E_max:.1f} kJ/mol  ({E_max/kBT:.1f} kT)")

# Estimate where P drops to 1/e, 1/e^2 etc
print(f"\nClassical wall estimate:")
for n_kt in [5, 10, 20]:
    r_wall = r0 + np.sqrt(2 * n_kt * kBT / kb)
    print(f"  G = {n_kt} kT  →  r = {r_wall*10:.4f} Å")

print(f"\nRecommended RETIS interfaces (reachable at {T:.0f}K):")
r_5kt  = r0 + np.sqrt(2 * 5  * kBT / kb)
r_10kt = r0 + np.sqrt(2 * 10 * kBT / kb)
print(f"  Safe range:    {r0*10:.3f} – {r_5kt*10:.4f} Å  (< 5 kT)")
print(f"  Pushable max:  up to ~{r_10kt*10:.4f} Å  (10 kT, rare but RETIS can reach)")
