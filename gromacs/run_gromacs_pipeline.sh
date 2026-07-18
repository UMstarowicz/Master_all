#!/bin/bash
# =============================================================================
# run_gromacs_pipeline.sh
# Full GROMACS pipeline: Energy Min → NVT equilibration → NPT production
#
# Usage:
#   bash run_gromacs_pipeline.sh spcfw
#   bash run_gromacs_pipeline.sh fbaeps
#   bash run_gromacs_pipeline.sh tip4p2005f
#
# Required files in the run directory before starting:
#   conf.gro      — starting structure
#   topol.top     — topology (with all #include paths correct)
#   em.mdp        — energy minimisation parameters
#   nvt.mdp       — NVT equilibration parameters
#   npt.mdp       — NPT production parameters
#
# Output:
#   em/           — energy minimisation results
#   nvt/          — NVT equilibration results
#   npt/          — NPT production trajectory + analysis inputs
# =============================================================================

set -euo pipefail   # exit on error, undefined var, pipe failure

MODEL=${1:-"model"}
GMX="gmx"
NTMPI=1
NTOMP=8
MAXWARN=2

# Detect GPU flags (use if available)
GPU_FLAGS="-nb gpu -bonded gpu -pme gpu"

echo "============================================================"
echo "  GROMACS NPT pipeline — model: ${MODEL}"
echo "  $(date)"
echo "============================================================"

# Check required input files
for f in conf.gro topol.top em.mdp nvt.mdp npt.mdp; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: required file '$f' not found in $(pwd)"
        exit 1
    fi
done

# --------------------------------------------------------------------------- #
#  STEP 1 — Energy minimisation                                               #
# --------------------------------------------------------------------------- #
echo ""
echo ">>> STEP 1: Energy minimisation"
mkdir -p em
cd em

${GMX} grompp -f ../em.mdp \
              -c ../conf.gro \
              -p ../topol.top \
              -o em.tpr \
              -maxwarn ${MAXWARN}

${GMX} mdrun -v \
             -ntmpi ${NTMPI} \
             -ntomp ${NTOMP} \
             -deffnm em

# Quick check: did EM converge?
FMAX=$(grep "Maximum force" em.log 2>/dev/null | tail -1 | awk '{print $NF}' || echo "unknown")
echo "  EM final Fmax: ${FMAX} kJ/mol/nm"

cd ..
echo "  EM done — minimised structure: em/em.gro"

# --------------------------------------------------------------------------- #
#  STEP 2 — NVT equilibration                                                 #
# --------------------------------------------------------------------------- #
echo ""
echo ">>> STEP 2: NVT equilibration (200 ps)"
mkdir -p nvt
cd nvt

${GMX} grompp -f ../nvt.mdp \
              -c ../em/em.gro \
              -p ../topol.top \
              -o nvt.tpr \
              -maxwarn ${MAXWARN}

${GMX} mdrun -v \
             -ntmpi ${NTMPI} \
             -ntomp ${NTOMP} \
             ${GPU_FLAGS} \
             -deffnm nvt || \
${GMX} mdrun -v \
             -ntmpi ${NTMPI} \
             -ntomp ${NTOMP} \
             -deffnm nvt   # retry without GPU flags if GPU unavailable

cd ..
echo "  NVT done — equilibrated structure: nvt/nvt.gro"

# --------------------------------------------------------------------------- #
#  STEP 3 — NPT production                                                    #
# --------------------------------------------------------------------------- #
echo ""
echo ">>> STEP 3: NPT production (10 ns)"
mkdir -p npt
cd npt

${GMX} grompp -f ../npt.mdp \
              -c ../nvt/nvt.gro \
              -t ../nvt/nvt.cpt \
              -p ../topol.top \
              -o npt.tpr \
              -maxwarn ${MAXWARN}

${GMX} mdrun -v \
             -ntmpi ${NTMPI} \
             -ntomp ${NTOMP} \
             ${GPU_FLAGS} \
             -deffnm npt || \
${GMX} mdrun -v \
             -ntmpi ${NTMPI} \
             -ntomp ${NTOMP} \
             -deffnm npt

cd ..
echo "  NPT done — trajectory: npt/npt.xtc"

# --------------------------------------------------------------------------- #
#  STEP 4 — Quick sanity checks on NPT run                                   #
# --------------------------------------------------------------------------- #
echo ""
echo ">>> STEP 4: Sanity checks"
cd npt

# Average density
echo "  Density (should be ~995–1000 kg/m³ at 300 K, 1 bar):"
echo "22 0" | ${GMX} energy -f npt.edr -o density_${MODEL}.xvg -quiet 2>/dev/null || \
echo "    (run 'gmx energy -f npt/npt.edr' manually to check density)"

# Average temperature
echo "  Temperature:"
echo "15 0" | ${GMX} energy -f npt.edr -o temperature_${MODEL}.xvg -quiet 2>/dev/null || \
echo "    (run 'gmx energy -f npt/npt.edr' manually to check temperature)"

# Average pressure
echo "  Pressure:"
echo "16 0" | ${GMX} energy -f npt.edr -o pressure_${MODEL}.xvg -quiet 2>/dev/null || \
echo "    (run 'gmx energy -f npt/npt.edr' manually to check pressure)"

cd ..

echo ""
echo "============================================================"
echo "  Pipeline complete: ${MODEL}"
echo "  $(date)"
echo "  Trajectory for analysis: npt/npt.xtc"
echo "  Run next: python analyse_trajectory.py --model ${MODEL}"
echo "============================================================"
