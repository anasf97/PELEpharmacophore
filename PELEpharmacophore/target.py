import os
import re
import numpy as np
import glob
import copy
from scipy.spatial import distance
from multiprocessing import Pool
import PELEpharmacophore.grid as gr
import PELEpharmacophore.helpers as hl

class Target():

    def __init__(self, indir):
        self.result_dir = f"{indir}/"
        self.trajectories = glob.glob(os.path.join(self.result_dir, "trajectory_*.pdb"))
        self.reports = glob.glob(os.path.join(self.result_dir, "report_*"))
        self.match_traj_and_report()
        self.chain = None

    def match_traj_and_report(self):
        self.trajectories.sort()
        self.reports.sort()
        self.traj_and_reports = zip(self.trajectories, self.reports)

    def set_ligand(self, chain, name, residue):
        self.chain = chain
        self.name = name
        self.residue = residue

    def set_features(self, features):
        self.features = features

    def set_grid(self, center, radius):
        self.grid = gr.Grid(center, radius)
        self.grid.generate_voxels()

    def get_structure(self, file):
        structure = hl.read_pdb(file)
        return structure

    def get_grid_atoms(self, model):
        dist = np.sqrt(3)*self.grid.radius #get dist from center to vertex
        featured_atoms = []
        atoms = []
        for feature, atomlist in self.features.items():
            for atom in atomlist:
                bio_atom = model[self.chain][(f"H_{self.name}", self.residue, " ")][atom]
                atoms.append(bio_atom)

                a = Atom(bio_atom)
                a.set_feature(feature)
                featured_atoms.append(a)

        atoms_near = hl.neighbor_search(atoms, self.grid.center, dist)
        grid_atoms = [a for a in atoms_near if hl.inside_grid(a, self.grid.v1, self.grid.v8)]
        featured_grid_atoms = [fa for fa in featured_atoms if fa.atom in grid_atoms ]
        return featured_grid_atoms

    def check_voxels(self, grid_atoms):
        model_grid = copy.deepcopy(self.grid)
        voxel_centers = np.array([v.center for v in model_grid.voxels])
        atom_coords = np.array([a.atom.get_coord() for a in grid_atoms])
        dist = distance.cdist(atom_coords, voxel_centers, 'euclidean')
        min = dist.argmin(axis=1) #get index of closest voxel to each atom
        voxel_dict = {}
        for i, atom in enumerate(grid_atoms):
            feature = atom.get_feature()
            origin = atom.get_origin()
            voxel = min[i]
            model_grid.voxels[voxel].count_feature(feature) #ir añadiendo la frecuencia de cada atomo en un diccionario
            model_grid.voxels[voxel].add_origin(feature, origin) #añadir tambien DICCIONARIO con los modelos de origen
        return model_grid

    def analyze_trajectory(self, traj_and_report):
        trajfile, report = traj_and_report
        trajectory = self.get_structure(trajfile)
        accepted_steps = hl.accepted_pele_steps(report)
        traj_grid = copy.deepcopy(self.grid)
        for step in accepted_steps:
            model = trajectory[0]
            grid_atoms = self.get_grid_atoms(model)
            model_grid = self.check_voxels(grid_atoms)
            traj_grid = self.merge_grids(traj_grid, model_grid)
        return traj_grid

    def merge_grids(self, grid, other_grid):
        if grid.is_empty():
            grid = copy.deepcopy(other_grid)
        else:
            for i, other_voxel in enumerate(other_grid.voxels):
                voxel = grid.voxels[i]
                if other_voxel.freq_dict:
                    for feature, other_freq in other_voxel.freq_dict.items():
                        hl.frequency_dict(voxel.freq_dict, feature, other_freq)
                        other_models = other_voxel.origin_dict[feature]
                        [hl.list_dict(voxel.origin_dict, feature, model) for model in other_models]
        return grid

    def set_frequency_filter(self, threshold):
        freq_dict_all = {}
        self.threshold_dict = {}
        for voxel in self.grid.voxels:
            if voxel.freq_dict:
                for feature, freq in voxel.freq_dict.items():
                    hl.list_dict(freq_dict_all, feature, freq)
        for element, freqlist in freq_dict_all.items():
            hist, bin_edges = np.histogram(freqlist)
            self.threshold_dict[element] = bin_edges[threshold]

    def save_pharmacophores(self, outdir="Pharmacophores"):
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
                            feature_origin = voxel.origin_dict[feature]
                            f.write(hl.format_line_pdb(voxel.center, feature, freq, feature_origin))


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
            match = re.search("trajectory_(\d)*", trajectory)
            trajectory = int(match.group(1))
            model = self.atom.get_full_id()[1]
            self.origin = (trajectory, model)
        return self.origin


if __name__ == "__main__":
    target = Target("/home/ana/test4")
    print(target.trajectories)
    target.set_ligand("L", "FRA", 900)
    #features={'HBD': ['NC1'], 'HBA': ['NB1', 'NC3', 'O2'], 'ALI': ['FD3', 'C1'], 'ARO': ['CA5', 'CD1']}
    features={'NEG': ['C2'], 'ALI': ['C1']}
    target.set_features(features)
    target.set_grid((2.173, 15.561, 28.257), 7)
    with Pool(5) as p: # add n_workers as arg
        traj_grids = p.map(target.analyze_trajectory, target.traj_and_reports)
        p.close()
        p.join()

    for traj_grid in traj_grids:
        target.grid = target.merge_grids(target.grid, traj_grid)

    target.set_frequency_filter(2)
    target.save_pharmacophores()

    for v in target.grid.voxels:
        print(v.freq_dict)

    # s = Snap("processed_1a9u.pdb")
    # s.set_grid((2.173, 15.561, 28.257), 6)
    # a = s.get_grid_atoms()
    # import matplotlib.pyplot as plt
    # from mpl_toolkits.mplot3d import Axes3D
    #
    # x, y, z = [], [], []
    # for atom in a:
    #     x.append(atom.get_coord()[0])
    #     y.append(atom.get_coord()[1])
    #     z.append(atom.get_coord()[2])
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(x, y, z, c='r', marker='o')
    #
    # ax.set_xlabel('X Label')
    # ax.set_ylabel('Y Label')
    # ax.set_zlabel('Z Label')
    #
    # plt.show()
