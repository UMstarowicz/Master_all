Simulation
----------
task = retis
steps = 20000
interfaces = [0.11018, 0.111, 0.112, 0.113, 0.1135]
zero_left = 0.109

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
subcycles = 10
gmx_format = gro
maxwarn = 2

TIS settings
------------
freq = 0.5
maxlength = 10000
aimless = True
allowmaxlength = False
zero_momentum = False
rescale_energy = False
sigma_v = -1
seed = 0
shooting_move = sh
relative_shoots = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3]

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
module = orderp_grotthuss_cvs.py

Output settings
---------------
pathensemble-file = 1
order-file = 1
energy-file = 1
trajectory-file = 50
