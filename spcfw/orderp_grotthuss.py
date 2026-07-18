# -*- coding: utf-8 -*-
"""
Order parameter for water autoionization — Moqadam et al. PNAS 2018.
λ = largest covalent O-H bond distance in the system (nm).

Simplified: no PBC needed for intramolecular O-H bonds (~0.10 nm).
"""
import numpy as np
from pyretis.orderparameter import OrderParameter
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class WaterAutoionization(OrderParameter):
    """Largest O-H bond distance in system (nm). Moqadam et al. 2018."""

    def __init__(self):
        super().__init__(description='Largest O-H distance (nm)')
        self.oh_cutoff = 0.14   # nm — max covalent O-H distance

    def calculate(self, system):
        """
        Returns [lambda] where lambda = max O-H distance in nm.
        SPC/Fw atom order per molecule: OW, HW1, HW2 (indices 3i, 3i+1, 3i+2)
        """
        pos   = system.particles.pos   # shape (n_atoms, 3), units nm
        n_mol = pos.shape[0] // 3
        max_oh = 0.0

        for i in range(n_mol):
            o  = pos[3 * i]
            h1 = pos[3 * i + 1]
            h2 = pos[3 * i + 2]

            r1 = float(np.sqrt(np.sum((h1 - o) ** 2)))
            r2 = float(np.sqrt(np.sum((h2 - o) ** 2)))

            # Only covalent bonds (exclude H-bonds ~0.18 nm)
            if r1 < self.oh_cutoff and r1 > max_oh:
                max_oh = r1
            if r2 < self.oh_cutoff and r2 > max_oh:
                max_oh = r2

        return [max_oh]
