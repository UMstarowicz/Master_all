Retis Hydrate Nucleation
========================

Simulation
----------
task = retis
steps = 1250
interfaces = [12, 20, 30, 40, 50, 60, 70, 80, 90, 100,
              110, 120, 130, 140, 150,  180,
              210, 400]
zero_left = 10

System
------
units = gromacs

Engine settings
---------------
class = gromacs2
gmx = gmx
mdrun = gmx mdrun -ntmpi 1 -ntomp 8
input_path = gromacs_input
timestep = 0.002
subcycles = 5000
gmx_format = gro
maxwarn = 1

TIS settings
------------
freq = 0.5
maxlength = 100000
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
class = MCG
module = op_cv.py

Output settings
---------------
pathensemble-file = 1
order-file = 1
energy-file = 1
trajectory-file = 100
