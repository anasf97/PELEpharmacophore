import os
import numpy as np
from scipy.spatial import distance
import glob
from multiprocessing import Pool
import PELEpharmacophore.grid as gr
import PELEpharmacophore.helpers as hl

class Target():

    def __init__(self, indir):
        self.filelist = glob.glob(f"{indir}*pdb")
        self.voxel_dict_all = {}
        self.chain = None

    def set_target_chain(self, target_chain):
        self.target_chain = target_chain

    def set_ligand(self, chain, name, residue):
        self.chain = chain
        self.name = name
        self.residue = residue

    def set_features(self, features):
        self.features = features

    def set_grid(self, center, radius=7):
        self.grid = gr.Grid(center, radius)
        self.grid.generate_voxels()

    def set_atom_ref(self, atom):
        self.atom_ref = atom

    def get_structure(self, file):
        structure = hl.read_pdb(file)
        return structure

    def get_side_chain_atoms_around_ligand(self, structure, dist=4):
        atoms = []
        for model in structure:
            model_atoms = []
            target_atoms = []
            [target_atoms.append(a) for a in model[self.target_chain].get_atoms() if a.get_name() not in ("CA","N","C")]
            ligand = model[self.chain][(f"H_{self.name}", self.residue, " ")]
            #atom_ref = model[self.chain][(f"H_{self.name}", self.residue, " ")][self.atom_ref]
            #center = atom_ref.get_coord()
            for atom in ligand.get_atoms():
                atoms_near = hl.neighbor_search(target_atoms, atom.get_coord(), dist)
                [model_atoms.append(a) for a in atoms_near]
            model_atoms = set(model_atoms)
        [atoms.append(a) for a in model_atoms]
        return atoms

    def get_grid_atoms(self, structure):
        dist = np.sqrt(3)*self.grid.radius #get dist from center to vertex
        featured_atoms = []
        atoms = []
        for model in structure:
            for feature, atomlist in self.features.items():
                for atom in atomlist:
                    bio_atom = model[self.chain][(f"H_{self.name}", self.residue, " ")][atom]
                    atoms.append(bio_atom)

                    a = FeaturedAtom(bio_atom)
                    a.set_feature(feature)
                    featured_atoms.append(a)

        atoms_near = hl.neighbor_search(atoms, self.grid.center, dist)
        grid_atoms = [a for a in atoms_near if hl.atoms_inside_grid(a, self.grid.v1, self.grid.v8)]
        featured_grid_atoms = [fa for fa in featured_atoms if fa.atom in grid_atoms ]
        return featured_grid_atoms

    def check_voxels(self, atoms):
        voxel_centers = np.array([v.center for v in self.grid.voxels])
        atom_coords = np.array([a.get_coord() for a in atoms ])
        dist = distance.cdist(atom_coords, voxel_centers, 'euclidean')
        min = dist.argmin(axis=1) #get index of closest voxel to each atom
        voxel_dict = {}
        for i, atom in enumerate(atoms):
            #feature = atom.get_feature()
            element = atom.get_name()[0]
            voxel = min[i]
            voxel_dict.setdefault(voxel, []).append(element)
        return voxel_dict

    def analyze_trajectory(self, file):
        model = self.get_structure(file)
        #grid_atoms = self.get_grid_atoms(model)
        atoms_near = self.get_side_chain_atoms_around_ligand(model)
        voxel_dict = self.check_voxels(atoms_near)
        return voxel_dict

    def merge_voxel_dicts(self, voxel_dict):
        for voxel, atomlist in voxel_dict.items():
            [self.voxel_dict_all.setdefault(voxel, []).append(a) for a in atomlist]

    def get_frequencies(self):
        for i, atomlist in self.voxel_dict_all.items():
            freq_dict = {}
            for atom in atomlist:
                if (atom in freq_dict):
                    freq_dict[atom] += 1
                else:
                    freq_dict[atom] = 1
            self.grid.voxels[i].set_frequencies(freq_dict)

    def set_frequency_filter(self, threshold):
        freq_dict_all = {}
        self.threshold_dict = {}
        for voxel in self.grid.voxels:
            for element, freq in voxel.freq_dict.items():
                freq_dict_all.setdefault(element, []).append(freq)
        for element, freqlist in freq_dict_all.items():
            hist, bin_edges = np.histogram(freqlist)
            self.threshold_dict[element] = bin_edges[threshold]

    def save_pharmacophores(self, outdir="PharmacophoresTarget"):
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        for feature in self.threshold_dict:
            f = open(f"{outdir}/{feature}pharmacophore.pdb", 'w')
            f.close()
        for voxel in self.grid.voxels:
            for feature, freq in voxel.freq_dict.items():
                if freq >= self.threshold_dict[feature]:
                    with open(f"{outdir}/{feature}pharmacophore.pdb", 'a') as f:
                        f.write(hl.format_line_pdb(voxel.center, feature, freq))


class FeaturedAtom:

    def __init__(self, atom):
        self.atom = atom

    def set_feature(self, feature):
        self.feature = feature

    def get_feature(self):
        return self.feature


if __name__ == "__main__":
    target = Target("/home/ana/GitRepositories/PELEpharmacophore/PELEpharmacophore/0/trajectory_1*")
    target.set_ligand("L", "SB2", 800)
    target.set_target_chain("A")
    #target.set_features(features)
    #target.set_atom_ref("CC2")
    target.set_grid((2.173, 15.561, 28.257), 10)
    p = Pool(50)
    dicts = p.map(target.analyze_trajectory, target.filelist)
    p.close()
    p.join()
    for d in dicts:
        target.merge_voxel_dicts(d)
    target.get_frequencies()
    target.set_frequency_filter(1)
    target.save_pharmacophores()

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
