import os
import abc
import re
import glob
import numpy as np
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
        self.result_dir = f"{indir}/output/"
        self.top_file = os.path.join(self.result_dir, "topologies", "topology_0.pdb")
        self.trajectories = glob.glob(os.path.join(self.result_dir, "0",  "trajectory_*.pdb"))
        print(self.trajectories)
        print(self.top_file)
        self.reports = glob.glob(os.path.join(self.result_dir, "0", "report_*"))
        self.match_traj_and_report()
        self.chain = None


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
             Keys define the features and values, atoms associated with said feature.

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



    def get_indices(self, topology, resname):
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
        indices_dict = {}

        for feature, atomlist in self.features.items():
            indices = hl.get_indices(topology, resname, atomlist)
            indices_dict[feature] = indices
        return indices_dict

    @staticmethod
    def get_coordinates(traj_and_report, indices_dict):
        trajfile, report = traj_and_report
        indices = np.concatenate([i for i in indices_dict.values()])
        accepted_steps = hl.accepted_pele_steps(report)

        traj = hl.load_trajectory(trajfile, indices)
        coords = traj.xyz *10  # coord units from nm to A
        coords = coords[accepted_steps] # duplicate rows when a step is rejected

        coord_dict = {}
        start = 0
        for feature, indices in indices_dict.items():
            stop = start + len(indices)
            coord_dict[feature] = coords[:, start:stop, :].reshape(-1, 3)
            start += len(indices)
        return coord_dict

    @abc.abstractmethod
    def run(self):
        pass


    @abc.abstractmethod
    def save_pharmacophores(self):
        pass


class Feature:

    def __init__(self, atoms, feature):
        self.atoms = atoms
        self.feature = feature

        self.retrieve_origin()

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, atoms):
        if isinstance(atoms, tuple):
            self._atoms = atoms
        else:
            self._atoms = tuple([atoms])

    @property
    def feature(self):
        return self._feature

    @feature.setter
    def feature(self, feature):
        self._feature = feature

    @property
    def origin(self):
        return self._origin

    @origin.setter
    def origin(self, trajectory_model):
        self._origin = trajectory_model

    def retrieve_origin(self):
        atoms = self.atoms
        trajectory = hl.basename_without_extension(atoms[0].get_full_id()[0])
        match = re.search("trajectory_(\d+)*", trajectory)
        trajectory = int(match.group(1))
        model = atoms[0].get_full_id()[1]
        self.origin = (trajectory, model)

    def coordinates(self):
        atoms = self.atoms
        if len(atoms) == 1:
            return atoms[0].get_coord()
        if len(atoms) == 2:
            point1, point2 = (a.get_coord() for a in atoms)
            center = hl.midpoint(point1, point2)
            return center
        if len(atoms) == 3:
            point1, point2, point3 = (a.get_coord() for a in atoms)
            point4 = hl.midpoint(point1, point2)
            center = hl.midpoint(point3, point4)
            return center
