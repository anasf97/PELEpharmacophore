import os
import abc
import re
import glob
import numpy as np
import mdtraj as md
from scipy.spatial import distance
from collections import OrderedDict
from itertools import accumulate, chain
import PELEpharmacophore.helpers as hl
import PELEpharmacophore.data.fragment_features as ff


class SimulationAnalyzer(metaclass=abc.ABCMeta):
    """
    Class for analysing PELE simulations.
    """

    def __init__(self, indir, features=None):
        """
        Create a new SimulationAnalyzer object.

        Parameters
        ----------
        indir : str
             Name of the simulation directory.
        """
        subdirs = [os.path.join(indir, subdir) for subdir in os.listdir(indir)]
        output_dir = os.path.join(indir, "output")

        if output_dir in subdirs:
            self.simulations = [Simulation(indir, features)]
        else:
            self.simulations = [Simulation(s, features) for s in subdirs]


    def set_ligand(self, chain, resname, resnum):
        """
        Set the parameters that define the ligand.

        Parameters
        ----------
        chain : str
             Ligand chain name.
        resname : str
             Ligand residue name.
        resnum : int
             Ligand residue number.
        """
        self.chain = chain
        self.resname = resname
        self.resnum = resnum


    def get_topology(self, file):
        """
        Parses a PDB file and returns a structure object.

        Parameters
        ----------
        file : str
            PDB file path.

        Returns
        ----------
        structure : Bio.PDB.Structure
            Biopython structure object.
        """
        return hl.load_topology(file)



    def get_indices(self, topology, resname, atomlist, first_index=None):
        """
        Gets all atoms defined in the `features` attribute.

        Parameters
        ----------
        model : Bio.PDB.Model
            Biopython model object.

        Returns
        ----------
        featured_grid_atoms : list of Atom objects
        """
        indlst  = [hl.get_indices(topology, resname, a) for a in atomlist]
        indices = np.concatenate([list(i) for i in indlst])
        #res_indices =[i-first_index for i in indices]
        lengths = np.array([len(i) for i in indlst])

        return (indices , lengths)

    def get_target_indices(self):
        simulation = self.simulations[0]
        trajectory = md.load(simulation.topfile)
        topology = trajectory.topology
        bondgraph = topology.to_bondgraph()

        n_indices = topology.select("protein and symbol N")
        n_coords = trajectory.xyz[:, n_indices, :]
        n_coords = n_coords.reshape(-1, 3)
        atoms_near, n_inds = hl.neighbor_search(n_coords, self.center, 20)
        n_indices = n_indices[n_inds]

        o_indices = topology.select("protein and symbol O")
        o_coords = trajectory.xyz[:, o_indices, :]
        o_coords = o_coords.reshape(-1, 3)
        atoms_near, o_inds = hl.neighbor_search(o_coords, self.center, 20)
        o_indices = o_indices[o_inds]

        o_and_n_ind = n_indices + o_indices
        h_indices = []
        atoms = [topology.atom(i) for i in o_and_n_ind]
        for atom in atoms:
            for neighbor in bondgraph.neighbors(atom):
                if neighbor.element == "hydrogen":
                    h_indices.append(neighbor.index)


        return (n_indices, o_indices, h_indices)


    def get_coords(self, ncpus, steps=None):

        final_coord_dict={}

        target_indices = self.get_target_indices()

        for simulation in self.simulations:

            topology = self.get_topology(simulation.topfile)

            res_indices = topology.select(f"resname {self.resname}")

            first_index = res_indices[0]

            indices_dict = OrderedDict([(feature, self.get_indices(topology, self.resname, atomlist)) for feature, atomlist in simulation.features.items()])

            coord_dicts = hl.parallelize(get_coordinates, simulation.traj_and_reports, ncpus, indices_dict=indices_dict, resname=self.resname, steps=steps, target_indices=target_indices)

            sim_coord_dict = hl.merge_array_dicts(*coord_dicts)

            final_coord_dict = hl.merge_array_dicts(final_coord_dict, sim_coord_dict)


        return final_coord_dict


    @abc.abstractmethod
    def save_pharmacophores(self):
        pass


def get_coordinates(traj_and_report, indices_dict, resname, steps=None, target_indices=None):
    trajfile, report = traj_and_report
    n_indices, o_indices, h_indices = target_indices

    indices = np.concatenate([i[0] for i in indices_dict.values()])
    temp = np.argsort(indices)
    order_indices = np.empty_like(temp)
    order_indices[temp] = np.arange(len(indices))
    accepted_steps = hl.accepted_pele_steps(report)

    coords = hl.get_coordinates_from_trajectory(resname, trajfile, indices_to_retrieve=indices)
    coords = coords[:, order_indices]
    coords = coords[accepted_steps] # duplicate rows when a step is rejected
    if steps is not None:
        coords = coords[:steps]

    coord_dict = {}
    start = 0
    for feature, (indices, lengths) in indices_dict.items():
        stop = start + len(indices)
        feature_coords = coords[:, start:stop, :]
        start = stop
        feature_coords = calc_cycle_centroids(feature_coords, lengths)

        if feature == "HBA":
            target_coords = hl.get_coordinates_from_trajectory(resname, trajfile, indices_to_retrieve=h_indices)
            dist = 2.3
        if feature == "HBD":
            n_and_o_indices = n_indices + o_indices
            target_coords = hl.get_coordinates_from_trajectory(resname, trajfile, indices_to_retrieve=n_and_o_indices)
            dist = 2.3
        if feature == "POS":
            target_coords = hl.get_coordinates_from_trajectory(resname, trajfile, indices_to_retrieve=o_indices)
            dist = 4
        if feature == "NEG":
            target_coords = hl.get_coordinates_from_trajectory(resname, trajfile, indices_to_retrieve=n_indices)
            dist = 4

        if feature in ("HBA", "HBD", "POS", "NEG"):
            feature_coords = calc_distances(feature_coords, target_coords, dist)

        coord_dict[feature] = feature_coords.reshape(-1,3)

    return coord_dict


def calc_cycle_centroids(coords, lengths):
    ind = np.where(lengths > 1)[0]
    if ind.size == 0:
        return coords

    final_coords = []
    acc = list(accumulate(lengths))

    for i, length in enumerate(lengths):
        start = 0 if i == 0 else acc[i-1]
        stop = acc[i]
        if length > 1:
            cycle_coords = coords[:, start:stop, None]
            c = hl.centroid(cycle_coords)
        else:
            c = coords[:, start, None]

        final_coords.append(c)

    final_coords = np.hstack(final_coords)

    return final_coords


def calc_distances(feature_coords, traj_coords, dist):
    dist_matrix = np.array([distance.cdist(feature_coords[i], traj_coords[i], 'sqeuclidean') for i in range(len(feature_coords))])
    mask = np.any((dist_matrix <= dist**2), axis=2)
    features_interacting = feature_coords[mask]
    return features_interacting


class Simulation():
    """docstring for Simulation."""

    def __init__(self, indir, features=None):
        if features is None:
            frag_regex = r".*(?P<frag>frag\d+).*$"
            frag = re.match(frag_regex, indir)['frag']
            self.features = ff.fragment_features[frag]

        else:
            self.features = features

        self.indir = indir
        self.output = f"{indir}/output/"
        self.topfile = os.path.join(self.output, "topologies", "topology_0.pdb")
        self.trajectories = glob.glob(os.path.join(self.output, "0",  "trajectory_*.pdb"))
        self.reports = glob.glob(os.path.join(self.output, "0", "report_*"))
        self.traj_and_reports = self.match_traj_and_report()

    def match_traj_and_report(self):
        """
        Match each trajectory with its respective report.
        """
        self.trajectories.sort()
        self.reports.sort()
        traj_and_reports = list(zip(self.trajectories, self.reports))
        return traj_and_reports

    def set_features(d, fragment_features=ff.fragment_features):
        print(type(d))
        frag_regex = r".*(?P<frag>frag\d+$)"
        frag = re.match(frag_regex, d)['frag']
        features = fragment_features[frag]
        return features
