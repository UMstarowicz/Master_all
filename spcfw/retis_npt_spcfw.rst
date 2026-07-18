Retis Water Autoionization - SPC/Fw  NPT 300 K 1 bar
======================================================

Simulation
----------
task = retis
steps = 20000
interfaces = [0.1085, 0.1090, 0.1095, 0.1100, 0.1105, 0.1108, 0.1110, 0.1112, 0.1115, 0.1117, 0.1120]
zero_left = 0.1065

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
relative_shoots = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2]

RETIS settings
--------------
swapfreq = 0.5
relative_shoots = None
nullmoves = True
swapsimul = True

Initial-path
------------
method = load
load_and_kick = True
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
