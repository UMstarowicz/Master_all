import numpy as np
from pyretis.orderparameter import OrderParameter
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

class WaterAutoionization(OrderParameter):
    """Largest O-H bond distance (nm) for SPC/Fw water."""
    def __init__(self):
        super().__init__(description='Max O-H distance (nm)')
        self.oh_cutoff = 0.17

    def calculate(self, system):
        pos   = system.particles.pos
        n_mol = pos.shape[0] // 3
        max_oh = 0.0
        for i in range(n_mol):
            o  = pos[3*i]
            h1 = pos[3*i + 1]
            h2 = pos[3*i + 2]
            r1 = float(np.linalg.norm(h1 - o))
            r2 = float(np.linalg.norm(h2 - o))
            if r1 < self.oh_cutoff and r1 > max_oh:
                max_oh = r1
            if r2 < self.oh_cutoff and r2 > max_oh:
                max_oh = r2
        return [max_oh]
