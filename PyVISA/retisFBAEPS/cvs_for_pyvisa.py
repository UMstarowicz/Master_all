import numpy as np
import logging
from collections import deque
from pyretis.orderparameter import OrderParameter

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class MoqadamCVs(OrderParameter):
    """Collective variables from Moqadam et al. (2018) PNAS 115:E4569.

    Returns [w4, q, na, nd, qcos] for the water autoionization system.

    All positions are in nm (GROMACS units, SPC/Fw water O H H ordering).

    w4   : length of the 4-molecule hydrogen bond wire (nm)
    q    : tetrahedral order parameter around O_lambda (0-1)
    na   : H-bonds accepted by O_lambda's water molecule
    nd   : H-bonds donated by O_lambda's water molecule
    qcos : smallest cosine of the two internal wire angles
    """

    OO_CUT  = 0.350
    OH_CUT  = 0.245
    ANG_CUT = 30.0
    OH_COV  = 0.17
    WIRE_N  = 4

    def __init__(self):
        super().__init__(description='Moqadam CVs: w4, q, na, nd, qcos')

    def _mic(self, a, b, box):
        d = a - b
        d -= box * np.round(d / box)
        return d

    def _dist(self, a, b, box):
        return float(np.linalg.norm(self._mic(a, b, box)))

    def _find_o_lambda(self, pos, n_mol, box):
        max_oh, o_lam = -1.0, 0
        for i in range(n_mol):
            for h in (3*i+1, 3*i+2):
                d = self._dist(pos[3*i], pos[h], box)
                if d < self.OH_COV and d > max_oh:
                    max_oh, o_lam = d, 3*i
        return o_lam

    def _hbonds(self, pos, n_mol, box):
        hb = []
        for i in range(n_mol):
            do = 3*i
            for h in (3*i+1, 3*i+2):
                for j in range(n_mol):
                    ao = 3*j
                    if ao == do:
                        continue
                    if self._dist(pos[do], pos[ao], box) > self.OO_CUT:
                        continue
                    if self._dist(pos[h],  pos[ao], box) > self.OH_CUT:
                        continue
                    v1 = self._mic(pos[h],  pos[do], box)
                    v2 = self._mic(pos[ao], pos[do], box)
                    cos = np.dot(v1, v2) / (
                        np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-15)
                    ang = np.degrees(np.arccos(np.clip(cos, -1, 1)))
                    if ang < self.ANG_CUT:
                        hb.append((do, h, ao))
        return hb

    def _graph(self, hbonds, n_mol):
        g = {3*i: set() for i in range(n_mol)}
        for (d, h, a) in hbonds:
            g[d].add(a)
            g[a].add(d)
        return g

    def _wire(self, graph, start):
        q = deque([[start]])
        seen = set()
        while q:
            path = q.popleft()
            if len(path) == self.WIRE_N:
                return path
            key = tuple(path)
            if key in seen:
                continue
            seen.add(key)
            for nb in graph.get(path[-1], []):
                if nb not in path:
                    q.append(path + [nb])
        return None

    def _tet_q(self, centre, neigh4, box):
        vecs = np.array([self._mic(n, centre, box) for n in neigh4])
        vhat = vecs / (np.linalg.norm(vecs, axis=1, keepdims=True) + 1e-15)
        s = 0.0
        for j in range(3):
            for k in range(j+1, 4):
                s += (np.dot(vhat[j], vhat[k]) + 1/3)**2
        return float(1.0 - (3/8)*s)

    def calculate(self, system):
        """Return [w4, q, na, nd, qcos]."""
        pos   = system.particles.pos
        box   = system.box.length
        n_mol = len(pos) // 3

        o_lam  = self._find_o_lambda(pos, n_mol, box)
        hbonds = self._hbonds(pos, n_mol, box)
        graph  = self._graph(hbonds, n_mol)

        na = sum(1 for (_, _, a) in hbonds if a == o_lam)
        nd = sum(1 for (d, _, _) in hbonds if d == o_lam)

        wire = self._wire(graph, o_lam)
        if wire is not None:
            wp  = [pos[i] for i in wire]
            w4  = sum(self._dist(wp[i], wp[i+1], box) for i in range(3))
            cos = []
            for c in (1, 2):
                v1 = self._mic(wp[c-1], wp[c], box)
                v2 = self._mic(wp[c+1], wp[c], box)
                cos.append(float(np.clip(
                    np.dot(v1, v2) / (
                        np.linalg.norm(v1)*np.linalg.norm(v2) + 1e-15),
                    -1, 1)))
            qcos = min(cos)
        else:
            w4, qcos = 0.0, 0.0

        o_idx = [3*i for i in range(n_mol) if 3*i != o_lam]
        dists = [self._dist(pos[o_lam], pos[oi], box) for oi in o_idx]
        n4    = np.array([pos[o_idx[i]] for i in np.argsort(dists)[:4]])
        q     = self._tet_q(pos[o_lam], n4, box)

        return [w4, q, float(na), float(nd), qcos]
