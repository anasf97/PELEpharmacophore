import os
import abc
import re
import glob
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
        self.result_dir = f"{indir}/output/0"
        self.trajectories = glob.glob(os.path.join(self.result_dir, "trajectory_*.pdb"))
        self.reports = glob.glob(os.path.join(self.result_dir, "report_*"))
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


    def get_structure(self, file):
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
        structure = hl.read_pdb(file)
        return structure


    def get_atoms(self, model):
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
        featured_atoms = []
        for feature, atomlist in self.features.items():
            for atoms in atomlist:
                if isinstance(atoms, tuple):
                    at = tuple(model[self.chain][(f"H_{self.resname}", self.resnum, " ")][a] for a in atoms)
                else:
                    at = model[self.chain][(f"H_{self.resname}", self.resnum, " ")][atoms]

                f = Feature(at, feature)
                featured_atoms.append(f)
        return featured_atoms


    @abc.abstractmethod
    def analyze_trajectory(self):
        pass


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
