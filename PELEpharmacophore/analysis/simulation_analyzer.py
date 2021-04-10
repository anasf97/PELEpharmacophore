import os
import abc
import re
import glob
import numpy as np
from itertools import accumulate
import PELEpharmacophore.helpers as hl

class SimulationAnalyzer(metaclass=abc.ABCMeta):
    """
    Class for analysing PELE simulations.
    """

    def __init__(self, indir=None):
        """
        Create a new SimulationAnalyzer object.

        Parameters
        ----------
        indir : str
             Name of the simulation directory.
        """
        self.set_dir(indir)

    def set_dir(self, indir):
        self.result_dir = f"{indir}/output/"
        self.top_file = os.path.join(self.result_dir, "topologies", "topology_0.pdb")
        self.trajectories = glob.glob(os.path.join(self.result_dir, "0",  "trajectory_*.pdb"))
        self.reports = glob.glob(os.path.join(self.result_dir, "0", "report_*"))
        self.match_traj_and_report()

    def match_traj_and_report(self):
        """
        Match each trajectory with its respective report.
        """
        self.trajectories.sort()
        self.reports.sort()
        self.traj_and_reports = list(zip(self.trajectories, self.reports))


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


    def set_features(self, features):
        """
        Set the pharmacophore features of the ligand.

        Parameters
        ----------
        features : dict
             Dictionary of ligand features.
             Keys define the features and values, the atoms associated with said feature.

        Examples
        ----------
        >>> features = {'HBD': ['NC1'], 'HBA': ['NB1', 'NC3', 'O2']}
        """
        self.features = features


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



    def get_indices(self, topology, resname, atomlist):
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
        lengths = np.array([len(i) for i in indlst])

        return (indices , lengths)


    @abc.abstractmethod
    def run(self, ncpus):
        import tracemalloc

        tracemalloc.start()

        topology = self.get_topology(self.top_file)

        indices_dict = {feature: self.get_indices(topology, self.resname, atomlist) \
                        for feature, atomlist in self.features.items()}
        first_size, first_peak = tracemalloc.get_traced_memory()

        coord_dicts = hl.parallelize(get_coordinates, self.traj_and_reports, ncpus, indices_dict=indices_dict)

        second_size, second_peak = tracemalloc.get_traced_memory()

        merged_coord_dict = hl.merge_array_dicts(*coord_dicts)


        print(f"First memory usage is {first_size / 10**6}MB; Peak was {first_peak / 10**6}MB")
        print(f"Second memory usage is {second_size / 10**6}MB; Peak was {second_peak / 10**6}MB")

        return merged_coord_dict


    @abc.abstractmethod
    def save_pharmacophores(self):
        pass

def get_coordinates(traj_and_report, indices_dict):
    trajfile, report = traj_and_report
    indices = np.concatenate([i[0] for i in indices_dict.values()])
    accepted_steps = hl.accepted_pele_steps(report)

    traj = hl.load_trajectory(trajfile, indices)
    coords = traj.xyz *10  # coord units from nm to A
    coords = coords[accepted_steps] # duplicate rows when a step is rejected

    coord_dict = {}
    start = 0
    for feature, (indices, lengths) in indices_dict.items():
        stop = start + len(indices)
        feature_coords = coords[:, start:stop, :]
        start = stop
        feature_coords = calc_cycle_centroids(feature_coords, lengths)
        coord_dict[feature] = feature_coords.reshape(-1, 3)
    return coord_dict

def calc_cycle_centroids(coords, lengths):
    ind = np.where(lengths > 1)[0]
    if ind.size == 0:
        return coords

    acc = list(accumulate(lengths))
    for i in ind:
        start = 0 if i == 0 else acc[i-1]
        stop = acc[i]
        cycle_coords = coords[:, start:stop, :]
        centroid = hl.centroid(cycle_coords)
        coords[:, start, :] = centroid
        coords[:, start+1:stop, :] = np.NaN
    coords = coords[~np.all(np.isnan(coords), axis=2)]
    return coords
