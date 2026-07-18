# -*- coding: utf-8 -*-
# Copyright (c) 2023, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""
This file includes an order parameter for hydrate nucleation.

"""
import mdtraj as md
import numpy as np
import networkx as nx

import logging
from itertools import combinations
from pyretis.orderparameter import OrderParameter
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


class MCG(OrderParameter):
    """MCG(OrderParameter).

    This class computes the the MCG order parameter.

    Attributes
    ----------
    name : string
        A human-readable name for the order parameter.

    """

    def __init__(self):
        """Set up the order parameter.

        Parameters
        ----------
        index : tuple of ints
            This is the indices of the atom we will use the position of.

        """
        super().__init__(description='Water molecule distance')
        self.top = md.load_topology("gromacs_input/conf.gro")

    def calculate(self, system):
        """Calculate the order parameter.

        Here, the order parameter is just the distance between two
        particles.

        Parameters
        ----------
        system : object like `System` from `pyretis.core.system`
            This object is used to import the file names required
            for mdtraj.

        Returns
        -------
        out : float
            The order parameter.

        """

        def largest_component(adj_matrix):
            g = nx.from_numpy_array(adj_matrix)  # Convert matrix to graph
            largest_cc = max(nx.connected_components(g), key=len)
            return len(largest_cc)

        snapshot = md.Trajectory(
                    system.particles.pos.reshape((1, self.top.n_atoms, 3)),
                    self.top
                    # unitcell_vectors=system.box.reshape((1, 3, 3))
                )

        waters = snapshot.topology.select('name O')
        guests = snapshot.topology.select('name CZ or name C2')
        traj = snapshot  # md.Trajectory(positions, topology)

        neighbour_list = []
        node_list = []
        for guest in guests:
            neighbour_list_in = []
            n_count = 0

            neighbours = md.compute_neighbors(traj, 0.9, [guest], guests)[0]
            for neighbour in neighbours:
                water_list = []
                water_neighbours = np.intersect1d(
                    md.compute_neighbors(traj, 0.6, [guest], waters)[0],
                    md.compute_neighbors(traj, 0.6, [neighbour], waters)[0],
                    assume_unique=True)
                for water in water_neighbours:
                    angle = np.degrees(
                            md.compute_angles(traj,
                                              [[guest, neighbour, water]])[0])
                    if angle[0] < 45:

                        water_list.append(water)
                if len(water_list) >= 5:
                    n_count += 1
                    neighbour_list_in.append(neighbour)
            if n_count >= 1:
                neighbour_list.append(neighbour_list_in)
                node_list.append(guest)
        if node_list:
            adjacency_matrix = np.zeros([len(node_list), len(node_list)])
            for i, (node_1, neighbours) in enumerate(zip(node_list, neighbour_list)):
                count = 0
                for j, node_2 in enumerate(node_list):
                    if node_2 in neighbours:
                        adjacency_matrix[i, j] = 1
                        count += 1
                        if count == len(neighbours):
                            break
                    else:
                        adjacency_matrix[i, j] = 0
            mcg = largest_component(adjacency_matrix)
        else:
            mcg = 0
        return mcg
