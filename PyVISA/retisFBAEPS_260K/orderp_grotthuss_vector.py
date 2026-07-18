import numpy as np
from pyretis.orderparameter import OrderParameter
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

class WaterAutoionization(OrderParameter):
    """Largest O-H bond distance (nm) for SPC/Fw water (O H H O H H ...)."""

    def __init__(self):
        super().__init__(description='Max O-H distance (nm)')

    def calculate(self, system):
        pos = system.particles.pos          # shape (3*n_mol, 3)
        n_mol = pos.shape[0] // 3

        # Reshape: (n_mol, 3_atoms, 3_coords)
        mol_pos = pos[:3 * n_mol].reshape(n_mol, 3, 3)

        o  = mol_pos[:, 0, :]   # oxygen positions
        h1 = mol_pos[:, 1, :]   # first hydrogen
        h2 = mol_pos[:, 2, :]   # second hydrogen

        r1 = np.linalg.norm(h1 - o, axis=1)   # (n_mol,)
        r2 = np.linalg.norm(h2 - o, axis=1)

        return [float(np.maximum(r1, r2).max())]
