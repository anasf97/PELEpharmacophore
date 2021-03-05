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
        atoms = []
        for feature, atomlist in self.features.items():
            for atom in atomlist:
                bio_atom = model[self.chain][(f"H_{self.resname}", self.resnum, " ")][atom]
                atoms.append(bio_atom)

                a = Atom(bio_atom)
                a.set_feature(feature)
                featured_atoms.append(a)
        return atoms, featured_atoms


    @abc.abstractmethod
    def analyze_trajectory(self):
        pass


    @abc.abstractmethod
    def run(self):
        pass


    @abc.abstractmethod
    def save_pharmacophores(self):
        pass


class Atom:

    def __init__(self, atom):
        self.atom = atom
        self.origin = None

    def set_feature(self, feature):
        self.feature = feature

    def get_feature(self):
        return self.feature

    def get_origin(self):
        if self.origin is None:
            trajectory = hl.basename_without_extension(self.atom.get_full_id()[0])
            match = re.search("trajectory_(\d+)*", trajectory)
            trajectory = int(match.group(1))
            model = self.atom.get_full_id()[1]
            self.origin = (trajectory, model)
        return self.origin
