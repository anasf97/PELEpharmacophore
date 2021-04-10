import os
import copy
import numpy as np
from scipy.spatial import distance
import PELEpharmacophore.helpers as hl
import PELEpharmacophore.analysis.grid as gr
import PELEpharmacophore.analysis.simulation_analyzer as sa

class GridAnalyzer(sa.SimulationAnalyzer):
    """
    Class for analysing PELE simulations using a grid.
    """

    def set_grid(self, center, radius):
        """
        Set the grid .

        Parameters
        ----------
        center : list
            xyz coordinates of the grid center.
        radius : int
            Length of half the side of the grid.
        """
        self.grid = gr.Grid(center, radius)
        self.grid.generate_voxels()


    def get_grid_atoms(self, coords):
        """
        Gets all the atoms inside the grid for a given model.

        Parameters
        ----------
        model : Bio.PDB.Model
            Biopython model object.

        Returns
        ----------
        featured_grid_atoms : list of Atom objects
            Atoms inside the grid.
        """
        dist = np.sqrt(3)*self.grid.radius #get dist from center to vertex
        atoms_near = hl.neighbor_search(coords, self.grid.center, dist)
        inside_grid_mask = np.logical_and(self.grid.v1 <= atoms_near, atoms_near <=  self.grid.v8).all(axis=1)
        grid_atoms = atoms_near[inside_grid_mask]
        return grid_atoms


    def check_voxels(self, coords, voxel_centers):
        """
        Check which atoms are inside which voxel of the grid.

        Parameters
        ----------
        grid_atoms : list of Atom objects
            Atoms inside the grid.

        Returns
        ----------
        model_grid : Grid object
            Grid object that contains the given atoms in the correspondent voxels.
        """
        dist = distance.cdist(coords, voxel_centers, 'sqeuclidean')
        voxel_inds = dist.argmin(axis=1) #get index of closest voxel to each atom
        return voxel_inds


    def fill_grid(self, feature, voxel_inds):
        for i in voxel_inds:
            self.grid.voxels[i].count_feature(feature)


    def run(self, ncpus):
        """
        Analyze the full simulation.

        Parameters
        ----------
        ncpus : int
            Number of processors.
        """
        coord_dict = super().run(ncpus)

        voxel_centers = np.array([v.center for v in self.grid.voxels])

        grid_atoms_dict = {feature: self.get_grid_atoms(coords) \
                           for feature, coords in coord_dict.items()}

        voxel_ind_dict  = {feature: self.check_voxels(coords, voxel_centers) \
                           for feature, coords in grid_atoms_dict.items()}

        for feature, inds in voxel_ind_dict.items():
            self.fill_grid(feature, inds)


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
        freq_dict_all = {}
        self.threshold_dict = {}
        for voxel in self.grid.voxels:
            if voxel.freq_dict:
                for feature, freq in voxel.freq_dict.items():
                    freq_dict_all = hl.list_dict(freq_dict_all, feature, freq)
        for element, freqlist in freq_dict_all.items():
            hist, bin_edges = np.histogram(freqlist)
            self.threshold_dict[element] = bin_edges[threshold]
        return self.threshold_dict

    def save_pharmacophores(self, outdir="Pharmacophores"):
        """
        Save pharmacophore files in PDB format.

        Parameters
        ----------
        outdir : str
            Directory with the results.
        """
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        for feature in self.threshold_dict:
            path = os.path.join(outdir, f"{feature}pharmacophore.pdb")
            f = open(path, 'w')
            f.close()
        for voxel in self.grid.voxels:
            if voxel.freq_dict:
                for feature, freq in voxel.freq_dict.items():
                    path = os.path.join(outdir, f"{feature}pharmacophore.pdb")
                    if freq >= self.threshold_dict[feature]:
                        with open(path, 'a') as f:
                            f.write(hl.format_line_pdb(voxel.center, feature, freq))

if __name__ == "__main__":
    target = GridAnalyzer("/home/ana/GitRepositories/PELEpharmacophore/tests/data/simulation_1")
    target.set_ligand("L", "SB2", 800)
    features =  {'HBD': ['NC1'], 'HBA': ['NB1', 'NC3', 'O2'], 'ALI': ['FD3', 'C1'], 'ARO': [('CA1', 'CA4'), ('CD1', 'CD4'), ('CC4', 'CC5', 'CC2')]}
    target.set_features(features)
    target.set_grid((2.173, 15.561, 28.257), 7)
    target.run(1)
    target.set_frequency_filter(0)
