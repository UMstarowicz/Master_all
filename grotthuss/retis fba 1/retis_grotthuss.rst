; PyRETIS input for water autoionization / Grotthuss initiation detection
; Based on: Moqadam et al., PNAS 2018
;
; ORDER PARAMETER: λ = largest O-H bond distance in system (nm)
; Paper uses AIMD with BLYP/DZVP, 32 waters, 300 K, NVE
; We use SPC/Fw, 1000 waters, 300 K, NPT
;
; INTERFACES:
;   Paper (AIMD, Angstrom): 1.07, 1.10, 1.13, 1.16, 1.19, 1.22 ...
;   Here (nm, SPC/Fw):      0.107, 0.109, 0.111, 0.113, 0.115, 0.117, 0.120
;
;   SPC/Fw O-H equilibrium = 0.1012 nm, fluctuates ±0.005 nm typically
;   Interfaces cover the tail of the distribution up to stretched bonds
;   The "reactive initiation" threshold from paper = 1.15 Å = 0.115 nm
;
; KEY RESULT FROM PAPER (Fig. 3, 4, 5):
;   When λ crosses 1.15-1.16 Å (0.115-0.116 nm), the trajectory outcome
;   is predicted by: w4, q (tetrahedral OP), na (H-bonds accepted), qcos
;   The ML classifier (grotthuss_ml_classifier.py) computes these CVs
;   from the accepted PyRETIS paths.
;
; IMPORTANT NOTE:
;   With classical SPC/Fw, bonds cannot break → true autoionization
;   cannot occur. The simulation captures initiation conditions only.
;   For full reactive dynamics, replace GROMACS with CP2K+BLYP or
;   ASE+MACE-MP-0 (see documentation).

Retis Water Autoionization - Moqadam 2018
==========================================

Simulation
----------
task = retis
steps = 500
interfaces = [0.1085, 0.1090, 0.1095, 0.1100, 0.1105, 0.1110, 0.1115]
zero_left = 0.1078

System
------
units = gromacs

Engine settings
---------------
class = gromacs2
gmx = gmx
mdrun = gmx mdrun -ntmpi 1 -ntomp 8
input_path = gromacs_input
timestep = 0.001
subcycles = 5000
gmx_format = gro
maxwarn = 1

TIS settings
------------
freq = 0.5
maxlength = 50000
aimless = True
allowmaxlength = False
zero_momentum = False
rescale_energy = False
sigma_v = -1
seed = 0
shooting_move = 'sh'

RETIS settings
--------------
swapfreq = 0.5
relative_shoots = None
nullmoves = True
swapsimul = True

Initial-path
------------
method = load
load_folder = load

Orderparameter
--------------
class = WaterAutoionization
module = orderp_grotthuss.py

Output settings
---------------
pathensemble-file = 1
order-file = 1
energy-file = 1
trajectory-file = 50
