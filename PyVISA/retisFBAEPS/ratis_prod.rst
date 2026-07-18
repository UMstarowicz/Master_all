Simulation
----------
task = retis
steps = 10000
interfaces = [0.1070, 0.1100, 0.1130, 0.1160, 0.1190, 0.1220,
              0.1250, 0.1280, 0.1310, 0.1340, 0.1390, 0.1430,
              0.1480, 0.1520, 0.1560, 0.1800, 0.2000]
zero_left = 0.1050

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
subcycles = 10
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
relative_shoots = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3]

RETIS settings
--------------
swapfreq = 0.5
nullmoves = True
swapsimul = True

Initial-path
------------
method = load_and_kick
load = load/
kick-from = previous

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
