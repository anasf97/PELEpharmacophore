import os
import numpy as np
from itertools import chain
from sklearn.cluster import MeanShift
import PELEpharmacophore.helpers as hl
import PELEpharmacophore.analysis.simulation_analyzer as sa

class MeanshiftAnalyzer(sa.SimulationAnalyzer):
    """
    Class for analysing PELE simulations using the meanshift algorithm.
    """

    def analyze_trajectory(self, traj_and_report):
        """
        Analyze a given trajectory file.
        The analysis consist in, for each of the models in the trajectory,
        getting the atoms defined in the `features` attribute.

        Parameters
        ----------
        traj_and_report : tuple
            Trajectory and its respective report.

        Returns
        ----------
        all_featured_atoms : list of Atom objects
        """
        trajfile, report = traj_and_report
        trajectory = self.get_structure(trajfile)
        accepted_steps = hl.accepted_pele_steps(report)
        all_featured_atoms = []
        for step in accepted_steps:
            model = trajectory[step]
            atoms, featured_atoms = self.get_atoms(model)
            [all_featured_atoms.append(fa) for fa in featured_atoms]
        return all_featured_atoms

    def run(self, ncpus):
        """
        Analyze the full simulation.

        Parameters
        ----------
        ncpus : int
            Number of processors.
        """
        featured_atoms = hl.parallelize(self.analyze_trajectory, self.traj_and_reports, ncpus)
        featured_atoms = list(chain.from_iterable(featured_atoms))

        estimator = MeanShift(bandwidth=1, n_jobs=ncpus, cluster_all=True)
        atom_dict ={}
        for atom in featured_atoms:
            atom_dict = hl.list_dict(atom_dict, atom.get_feature(), atom.atom.get_coord())

        self.cluster_dict={}
        for feature, coords in atom_dict.items():
            results = estimator.fit_predict(coords)
            p_dict = {}
            for cluster in results:
                p_dict = hl.frequency_dict(p_dict, cluster, 1)

            for cluster, frequency in p_dict.items():
                center = estimator.cluster_centers_[cluster]
                c = Cluster(cluster, frequency, center)
                self.cluster_dict = hl.list_dict(self.cluster_dict, feature, c)

        return self.cluster_dict


    def set_frequency_filter(self, threshold):
        """
        Set a threshold to not show voxels with lower frequencies.
        For each feature, it creates a histogram of 10 bins with the frequencies of each voxel.

        Parameters
        ----------
        threshold : int
            Number of the bins that won't be saved.
            E. g.: if it's 1, the voxels on the first bin of the histogram won't be saved.
        """
        self.threshold_dict = {}
        for feature, clusters in self.cluster_dict.items():
            freqlist = [c.frequency for c in clusters]
            print(feature, freqlist)
            hist, bin_edges = np.histogram(freqlist)
            self.threshold_dict[feature] = bin_edges[threshold]
        print(self.threshold_dict)
        return self.threshold_dict

    def save_pharmacophores(self, outdir="Pharmacophores_ms"):
        """
        Save pharmacophore files in PDB format.

        Parameters
        ----------
        outdir : str
            Directory with the results.
        """
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        for feature in self.cluster_dict:
            path = os.path.join(outdir, f"{feature}pharmacophore.pdb")
            f = open(path, 'w')
            f.close()
        for feature, clusters in self.cluster_dict.items():
            path = os.path.join(outdir, f"{feature}pharmacophore.pdb")
            with open(path, 'a') as f:
                for cluster in clusters:
                    if cluster.frequency >= self.threshold_dict[feature]:
                        f.write(hl.format_line_pdb(cluster.center, feature, cluster.frequency))


class Cluster(object):
    """docstring for Cluster."""

    def __init__(self, id, frequency, center):
        self.id = id
        self.frequency = frequency
        self.center = center

if __name__ == "__main__":
    target = MeanshiftAnalyzer("/gpfs/scratch/bsc72/bsc72801/ana_sanchez/test5")
    target.set_ligand("L", "FRA", 900)
    #features={'HBD': ['NC1'], 'HBA': ['NB1', 'NC3', 'O2'], 'ALI': ['FD3', 'C1'], 'ARO': ['CA5', 'CD1']}
    features={'NEG': ['C2'], 'ALI': ['C1']}
    target.set_features(features)
    target.run(20)
    target.set_frequency_filter(0)
    target.save_pharmacophores("PharmacophoresTest5_ms")
