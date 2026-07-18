import numpy as np
from pyretis.orderparameter import OrderParameter
import logging

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class WaterAutoionization(OrderParameter):
    """
    CV set for proton transfer / Grotthuss analysis.

    Returns:
    OP1 = max O-H bond distance
    OP2 = wire compression w4
    OP3 = tetrahedrality q
    OP4 = accepted H-bonds na
    OP5 = wire alignment qcos
    """

    def __init__(self):
        super().__init__(
            description='''
            OP1: max O-H distance
            OP2: wire compression w4
            OP3: tetrahedrality q
            OP4: accepted H-bonds na
            OP5: wire alignment qcos
            '''
        )

        self.oh_cutoff = 0.17
        self.hbond_cutoff = 0.35

    def calculate(self, system):

        pos = system.particles.pos
        n_mol = pos.shape[0] // 3

        # ============================================================
        # OP1: largest O-H distance
        # ============================================================

        max_oh = 0.0

        oxygen_positions = []

        for i in range(n_mol):

            o = pos[3 * i]
            h1 = pos[3 * i + 1]
            h2 = pos[3 * i + 2]

            oxygen_positions.append(o)

            r1 = float(np.linalg.norm(h1 - o))
            r2 = float(np.linalg.norm(h2 - o))

            if r1 < self.oh_cutoff and r1 > max_oh:
                max_oh = r1

            if r2 < self.oh_cutoff and r2 > max_oh:
                max_oh = r2

        oxygen_positions = np.array(oxygen_positions)

        # ============================================================
        # OP2: wire compression w4
        # Simple proxy:
        # average nearest-neighbor O-O distance
        # ============================================================

        oo_distances = []

        for i in range(n_mol):

            oi = oxygen_positions[i]

            min_dist = 999.0

            for j in range(n_mol):

                if i == j:
                    continue

                oj = oxygen_positions[j]

                d = np.linalg.norm(oi - oj)

                if d < min_dist:
                    min_dist = d

            oo_distances.append(min_dist)

        w4 = float(np.mean(oo_distances))

        # ============================================================
        # OP3: tetrahedrality q
        # simplified local tetrahedral order
        # ============================================================

        q_values = []

        for i in range(n_mol):

            oi = oxygen_positions[i]

            neighbor_vectors = []

            for j in range(n_mol):

                if i == j:
                    continue

                oj = oxygen_positions[j]

                vec = oj - oi
                d = np.linalg.norm(vec)

                if d < 0.35:
                    neighbor_vectors.append(vec)

            if len(neighbor_vectors) >= 4:

                q_local = 0.0
                count = 0

                for a in range(4):
                    for b in range(a + 1, 4):

                        va = neighbor_vectors[a]
                        vb = neighbor_vectors[b]

                        cosa = np.dot(va, vb) / (
                            np.linalg.norm(va)
                            * np.linalg.norm(vb)
                        )

                        q_local += (cosa + 1.0 / 3.0) ** 2
                        count += 1

                q_local = 1.0 - (3.0 / 8.0) * q_local

                q_values.append(q_local)

        tetra_q = float(np.mean(q_values)) if q_values else 0.0

        # ============================================================
        # OP4: accepted H-bonds na
        # crude estimate from O-O cutoff
        # ============================================================

        accepted_hbonds = 0

        for i in range(n_mol):

            oi = oxygen_positions[i]

            for j in range(n_mol):

                if i == j:
                    continue

                oj = oxygen_positions[j]

                d = np.linalg.norm(oi - oj)

                if d < self.hbond_cutoff:
                    accepted_hbonds += 1

        accepted_hbonds = accepted_hbonds / n_mol

        # ============================================================
        # OP5: wire alignment qcos
        # proxy alignment metric
        # ============================================================

        alignment_values = []

        for i in range(n_mol - 1):

            vec = oxygen_positions[i + 1] - oxygen_positions[i]

            norm = np.linalg.norm(vec)

            if norm > 1e-8:

                qcos = vec[2] / norm

                alignment_values.append(qcos)

        wire_alignment = (
            float(np.mean(alignment_values))
            if alignment_values
            else 0.0
        )

        return [
            max_oh,
            w4,
            tetra_q,
            accepted_hbonds,
            wire_alignment
        ]
