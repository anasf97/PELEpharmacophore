import os
import numpy as np
from itertools import chain
from sklearn.cluster import MeanShift
import PELEpharmacophore.helpers as hl
import PELEpharmacophore.analysis.simulation_analyzer as sa
import PELEpharmacophore.output.pharmacophore_writer as pw

class MeanshiftAnalyzer(sa.SimulationAnalyzer):
    """
    Class for analysing PELE simulations using the meanshift algorithm.
    """

    def run(self, ncpus, steps=None):
        """
        Analyze the full simulation.

        Parameters
        ----------
        ncpus : int
            Number of processors.
        """
        coord_dict = self.get_coords(ncpus, steps)

        estimator = MeanShift(bandwidth=1, n_jobs=ncpus, cluster_all=True)

        self.cluster_dict={}
        for feature, coords in coord_dict.items():
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
            hist, bin_edges = np.histogram(freqlist)
            self.threshold_dict[feature] = bin_edges[threshold]

        for feature, clusters in self.cluster_dict.items():
            for i, cluster in enumerate(clusters):
                if cluster.frequency < self.threshold_dict[feature]:
                    clusters.pop(i)
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

        name = self.simulations[0].indir
        pw.PharmacophoreWriter(name, self.voxel_dict, self.coords, outdir)


class Cluster(object):
    """docstring for Cluster."""

    def __init__(self, id, frequency, center):
        self.id = id
        self.frequency = frequency
        self.center = center

if __name__ == "__main__":
    target = MeanshiftAnalyzer("/home/ana/GitRepositories/PELEpharmacophore/tests/data/simulation_1")
    target.set_ligand("L", "SB2", 800)
    #features={'HBD': ['NC1'], 'HBA': ['NB1', 'NC3', 'O2'], 'ALI': ['FD3', 'C1'], 'ARO': ['CA5', 'CD1']}
    #features={'NEG': ['C2'], 'ALI': ['C1']}
    features =  {'HBD': ['NC1'], 'HBA': ['NB1', 'NC3', 'O2'], 'ALI': ['FD3', 'C1'], 'ARO': [('CA1', 'CA4'), ('CD1', 'CD4'), ('CC4', 'CC5', 'CC2')]}
    target.set_features(features)
    target.run(1)
    for feature, clusters in target.cluster_dict.items():
        centers = [cluster.center for cluster in clusters]
        freqs = [cluster.frequency for cluster in clusters]

        print(feature, freqs)
    # target.set_frequency_filter(0)
    # target.save_pharmacophores("PharmacophoresTest5_ms")
