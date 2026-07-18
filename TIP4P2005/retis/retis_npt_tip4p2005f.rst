Retis Water Autoionization - TIP4P/2005f  NPT 300 K 1 bar
===========================================================
; Fixed from original retis_prod.rst:
;   - subcycles : 10    → 5000  (was 0.01 ps/segment — too short to sample)
;   - relative_shoots: 17 entries → 5 entries (must match n_interfaces)
;   - ref-t: 270 K → 300 K  (in MDP)
;   - pcoupl: no → C-rescale  (in MDP, for NPT)

Simulation
----------
task = retis
steps = 20000
interfaces = [0.107, 0.108, 0.109, 0.110, 0.1105]
zero_left = 0.105

System
------
units = gromacs

Engine settings
---------------
class = gromacs2
gmx = gmx
mdrun = gmx mdrun -ntmpi 1 -ntomp 8 -nb gpu -bonded gpu
input_path = gromacs_input
timestep = 0.001
subcycles = 5000
gmx_format = gro
maxwarn = 2

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
shooting_move = sh
relative_shoots = [0, 0, 0, 1, 1]

RETIS settings
--------------
swapfreq = 0.5
nullmoves = True
swapsimul = True

Initial-path
------------
method = load
load_and_kick = True
load_folder = load/

Orderparameter
--------------
class = WaterAutoionization
module = orderp_tip4p.py

Output settings
---------------
pathensemble-file = 1
order-file = 1
energy-file = 1
trajectory-file = 50
