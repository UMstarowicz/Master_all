PyRETIS input settings
======================
For more info, please see: http://www.pyretis.org
Have Fun!

Simulation settings
-------------------
task = 'retis'
steps = 20000
interfaces = [0.107, 0.108, 0.109, 0.11, 0.1105]
zero_left = 0.105
exe_path = '/home/ula/pyretis_run/project_masters/retis_tip4p2005f'
rgen = 'rgen'
priority_shooting = False
flux = True
zero_ensemble = True
endcycle = 20000

System settings
---------------
units = 'gromacs'
dimensions = 3
temperature = 1.0

Engine settings
---------------
class = 'gromacs2'
gmx = 'gmx'
mdrun = 'gmx mdrun -ntmpi 1 -ntomp 8 -nb gpu -bonded gpu'
input_path = 'gromacs_input'
timestep = 0.001
subcycles = 10
gmx_format = 'gro'
maxwarn = 2
exe_path = '/home/ula/pyretis_run/project_masters/retis_tip4p2005f'
rgen = 'rgen'
type = 'external'
input_files = {'conf': '/home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/conf.gro',
               'index': '/home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/index.ndx',
               'input_o': '/home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/grompp.mdp',
               'topology': '/home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/topol.top'}

Particles settings
------------------
type = 'external'

Orderparameter settings
-----------------------
class = 'WaterAutoionization'
module = 'orderp_tip4p.py'
name = 'Order Parameter'

Output settings
---------------
pathensemble-file = 1
order-file = 1
energy-file = 1
trajectory-file = 50
backup = 'append'
cross-file = 1
restart-file = 1
screen = 10

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
shooting_move = 'sh'
relative_shoots = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3]
high_accept = False
rgen = 'rgen'
shooting_moves = []
mirror_freq = 0
target_freq = 0
target_indices = []

Initial-path settings
---------------------
method = 'load'
load_and_kick = True
load_folder = 'load/'

RETIS settings
--------------
swapfreq = 0.5
nullmoves = True
swapsimul = True

Ensemble
--------
simulation task = retis
simulation steps = 20000
simulation interfaces = [0.105, 0.107, 0.107]
simulation zero_left = 0.105
simulation exe_path = /home/ula/pyretis_run/project_masters/retis_tip4p2005f
simulation rgen = rgen
simulation priority_shooting = False
simulation flux = True
simulation zero_ensemble = True
system units = gromacs
system dimensions = 3
system temperature = 1.0
engine class = gromacs2
engine gmx = gmx
engine mdrun = gmx mdrun -ntmpi 1 -ntomp 8 -nb gpu -bonded gpu
engine input_path = gromacs_input
engine timestep = 0.001
engine subcycles = 10
engine gmx_format = gro
engine maxwarn = 2
engine exe_path = /home/ula/pyretis_run/project_masters/retis_tip4p2005f
engine rgen = rgen
engine type = external
engine input_files conf = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/conf.gro
engine input_files input_o = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/grompp.mdp
engine input_files topology = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/topol.top
engine input_files index = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/index.ndx
tis freq = 0.5
tis maxlength = 10000
tis aimless = True
tis allowmaxlength = False
tis zero_momentum = False
tis rescale_energy = False
tis sigma_v = -1
tis seed = 0
tis shooting_move = sh
tis relative_shoots = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3]
tis high_accept = False
tis rgen = rgen
tis shooting_moves = []
tis mirror_freq = 0
tis target_freq = 0
tis target_indices = []
tis ensemble_number = 0
tis detect = 0.108
retis swapfreq = 0.5
retis nullmoves = True
retis swapsimul = True
initial-path method = load
initial-path load_and_kick = True
initial-path load_folder = load/
orderparameter class = WaterAutoionization
orderparameter module = orderp_tip4p.py
orderparameter name = Order Parameter
output pathensemble-file = 1
output order-file = 1
output energy-file = 1
output trajectory-file = 50
output backup = append
output cross-file = 1
output restart-file = 1
output screen = 10
heading = {'text': 'PyRETIS input settings\n======================\nFor more info, please see: http://www.pyretis.org\nHave Fun!'}
particles type = external
analysis blockskip = 1
analysis bins = 100
analysis maxblock = 1000
analysis maxordermsd = -1
analysis ngrid = 1001
analysis plot plotter = mpl
analysis plot output = png
analysis plot style = pyretis
analysis report = ['latex', 'rst', 'html']
analysis skipcross = 1000
analysis txt-output = txt.gz
analysis tau_ref_bin = []
analysis skip = 0
interface = 0.107

Ensemble
--------
simulation task = retis
simulation steps = 20000
simulation interfaces = [0.107, 0.107, 0.1105]
simulation zero_left = 0.105
simulation exe_path = /home/ula/pyretis_run/project_masters/retis_tip4p2005f
simulation rgen = rgen
simulation priority_shooting = False
simulation flux = True
simulation zero_ensemble = True
system units = gromacs
system dimensions = 3
system temperature = 1.0
engine class = gromacs2
engine gmx = gmx
engine mdrun = gmx mdrun -ntmpi 1 -ntomp 8 -nb gpu -bonded gpu
engine input_path = gromacs_input
engine timestep = 0.001
engine subcycles = 10
engine gmx_format = gro
engine maxwarn = 2
engine exe_path = /home/ula/pyretis_run/project_masters/retis_tip4p2005f
engine rgen = rgen
engine type = external
engine input_files conf = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/conf.gro
engine input_files input_o = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/grompp.mdp
engine input_files topology = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/topol.top
engine input_files index = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/index.ndx
tis freq = 0.5
tis maxlength = 10000
tis aimless = True
tis allowmaxlength = False
tis zero_momentum = False
tis rescale_energy = False
tis sigma_v = -1
tis seed = 0
tis shooting_move = sh
tis relative_shoots = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3]
tis high_accept = False
tis rgen = rgen
tis shooting_moves = []
tis mirror_freq = 0
tis target_freq = 0
tis target_indices = []
tis ensemble_number = 1
tis detect = 0.108
retis swapfreq = 0.5
retis nullmoves = True
retis swapsimul = True
initial-path method = load
initial-path load_and_kick = True
initial-path load_folder = load/
orderparameter class = WaterAutoionization
orderparameter module = orderp_tip4p.py
orderparameter name = Order Parameter
output pathensemble-file = 1
output order-file = 1
output energy-file = 1
output trajectory-file = 50
output backup = append
output cross-file = 1
output restart-file = 1
output screen = 10
heading = {'text': 'PyRETIS input settings\n======================\nFor more info, please see: http://www.pyretis.org\nHave Fun!'}
particles type = external
analysis blockskip = 1
analysis bins = 100
analysis maxblock = 1000
analysis maxordermsd = -1
analysis ngrid = 1001
analysis plot plotter = mpl
analysis plot output = png
analysis plot style = pyretis
analysis report = ['latex', 'rst', 'html']
analysis skipcross = 1000
analysis txt-output = txt.gz
analysis tau_ref_bin = []
analysis skip = 0
interface = 0.107

Ensemble
--------
simulation task = retis
simulation steps = 20000
simulation interfaces = [0.107, 0.108, 0.1105]
simulation zero_left = 0.105
simulation exe_path = /home/ula/pyretis_run/project_masters/retis_tip4p2005f
simulation rgen = rgen
simulation priority_shooting = False
simulation flux = True
simulation zero_ensemble = True
system units = gromacs
system dimensions = 3
system temperature = 1.0
engine class = gromacs2
engine gmx = gmx
engine mdrun = gmx mdrun -ntmpi 1 -ntomp 8 -nb gpu -bonded gpu
engine input_path = gromacs_input
engine timestep = 0.001
engine subcycles = 10
engine gmx_format = gro
engine maxwarn = 2
engine exe_path = /home/ula/pyretis_run/project_masters/retis_tip4p2005f
engine rgen = rgen
engine type = external
engine input_files conf = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/conf.gro
engine input_files input_o = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/grompp.mdp
engine input_files topology = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/topol.top
engine input_files index = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/index.ndx
tis freq = 0.5
tis maxlength = 10000
tis aimless = True
tis allowmaxlength = False
tis zero_momentum = False
tis rescale_energy = False
tis sigma_v = -1
tis seed = 0
tis shooting_move = sh
tis relative_shoots = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3]
tis high_accept = False
tis rgen = rgen
tis shooting_moves = []
tis mirror_freq = 0
tis target_freq = 0
tis target_indices = []
tis ensemble_number = 2
tis detect = 0.109
retis swapfreq = 0.5
retis nullmoves = True
retis swapsimul = True
initial-path method = load
initial-path load_and_kick = True
initial-path load_folder = load/
orderparameter class = WaterAutoionization
orderparameter module = orderp_tip4p.py
orderparameter name = Order Parameter
output pathensemble-file = 1
output order-file = 1
output energy-file = 1
output trajectory-file = 50
output backup = append
output cross-file = 1
output restart-file = 1
output screen = 10
heading = {'text': 'PyRETIS input settings\n======================\nFor more info, please see: http://www.pyretis.org\nHave Fun!'}
particles type = external
analysis blockskip = 1
analysis bins = 100
analysis maxblock = 1000
analysis maxordermsd = -1
analysis ngrid = 1001
analysis plot plotter = mpl
analysis plot output = png
analysis plot style = pyretis
analysis report = ['latex', 'rst', 'html']
analysis skipcross = 1000
analysis txt-output = txt.gz
analysis tau_ref_bin = []
analysis skip = 0
interface = 0.108

Ensemble
--------
simulation task = retis
simulation steps = 20000
simulation interfaces = [0.107, 0.109, 0.1105]
simulation zero_left = 0.105
simulation exe_path = /home/ula/pyretis_run/project_masters/retis_tip4p2005f
simulation rgen = rgen
simulation priority_shooting = False
simulation flux = True
simulation zero_ensemble = True
system units = gromacs
system dimensions = 3
system temperature = 1.0
engine class = gromacs2
engine gmx = gmx
engine mdrun = gmx mdrun -ntmpi 1 -ntomp 8 -nb gpu -bonded gpu
engine input_path = gromacs_input
engine timestep = 0.001
engine subcycles = 10
engine gmx_format = gro
engine maxwarn = 2
engine exe_path = /home/ula/pyretis_run/project_masters/retis_tip4p2005f
engine rgen = rgen
engine type = external
engine input_files conf = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/conf.gro
engine input_files input_o = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/grompp.mdp
engine input_files topology = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/topol.top
engine input_files index = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/index.ndx
tis freq = 0.5
tis maxlength = 10000
tis aimless = True
tis allowmaxlength = False
tis zero_momentum = False
tis rescale_energy = False
tis sigma_v = -1
tis seed = 0
tis shooting_move = sh
tis relative_shoots = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3]
tis high_accept = False
tis rgen = rgen
tis shooting_moves = []
tis mirror_freq = 0
tis target_freq = 0
tis target_indices = []
tis ensemble_number = 3
tis detect = 0.11
retis swapfreq = 0.5
retis nullmoves = True
retis swapsimul = True
initial-path method = load
initial-path load_and_kick = True
initial-path load_folder = load/
orderparameter class = WaterAutoionization
orderparameter module = orderp_tip4p.py
orderparameter name = Order Parameter
output pathensemble-file = 1
output order-file = 1
output energy-file = 1
output trajectory-file = 50
output backup = append
output cross-file = 1
output restart-file = 1
output screen = 10
heading = {'text': 'PyRETIS input settings\n======================\nFor more info, please see: http://www.pyretis.org\nHave Fun!'}
particles type = external
analysis blockskip = 1
analysis bins = 100
analysis maxblock = 1000
analysis maxordermsd = -1
analysis ngrid = 1001
analysis plot plotter = mpl
analysis plot output = png
analysis plot style = pyretis
analysis report = ['latex', 'rst', 'html']
analysis skipcross = 1000
analysis txt-output = txt.gz
analysis tau_ref_bin = []
analysis skip = 0
interface = 0.109

Ensemble
--------
simulation task = retis
simulation steps = 20000
simulation interfaces = [0.107, 0.11, 0.1105]
simulation zero_left = 0.105
simulation exe_path = /home/ula/pyretis_run/project_masters/retis_tip4p2005f
simulation rgen = rgen
simulation priority_shooting = False
simulation flux = True
simulation zero_ensemble = True
system units = gromacs
system dimensions = 3
system temperature = 1.0
engine class = gromacs2
engine gmx = gmx
engine mdrun = gmx mdrun -ntmpi 1 -ntomp 8 -nb gpu -bonded gpu
engine input_path = gromacs_input
engine timestep = 0.001
engine subcycles = 10
engine gmx_format = gro
engine maxwarn = 2
engine exe_path = /home/ula/pyretis_run/project_masters/retis_tip4p2005f
engine rgen = rgen
engine type = external
engine input_files conf = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/conf.gro
engine input_files input_o = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/grompp.mdp
engine input_files topology = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/topol.top
engine input_files index = /home/ula/pyretis_run/project_masters/retis_tip4p2005f/gromacs_input/index.ndx
tis freq = 0.5
tis maxlength = 10000
tis aimless = True
tis allowmaxlength = False
tis zero_momentum = False
tis rescale_energy = False
tis sigma_v = -1
tis seed = 0
tis shooting_move = sh
tis relative_shoots = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3]
tis high_accept = False
tis rgen = rgen
tis shooting_moves = []
tis mirror_freq = 0
tis target_freq = 0
tis target_indices = []
tis ensemble_number = 4
tis detect = 0.1105
retis swapfreq = 0.5
retis nullmoves = True
retis swapsimul = True
initial-path method = load
initial-path load_and_kick = True
initial-path load_folder = load/
orderparameter class = WaterAutoionization
orderparameter module = orderp_tip4p.py
orderparameter name = Order Parameter
output pathensemble-file = 1
output order-file = 1
output energy-file = 1
output trajectory-file = 50
output backup = append
output cross-file = 1
output restart-file = 1
output screen = 10
heading = {'text': 'PyRETIS input settings\n======================\nFor more info, please see: http://www.pyretis.org\nHave Fun!'}
particles type = external
analysis blockskip = 1
analysis bins = 100
analysis maxblock = 1000
analysis maxordermsd = -1
analysis ngrid = 1001
analysis plot plotter = mpl
analysis plot output = png
analysis plot style = pyretis
analysis report = ['latex', 'rst', 'html']
analysis skipcross = 1000
analysis txt-output = txt.gz
analysis tau_ref_bin = []
analysis skip = 0
interface = 0.11

Analysis settings
-----------------
blockskip = 1
bins = 100
maxblock = 1000
maxordermsd = -1
ngrid = 1001
plot = {'output': 'png', 'plotter': 'mpl', 'style': 'pyretis'}
report = ['latex', 'rst', 'html']
skipcross = 1000
txt-output = 'txt.gz'
tau_ref_bin = []
skip = 0