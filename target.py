from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch
import numpy as np
from scipy.spatial import distance
import grid as gr
import glob
from multiprocessing import Pool

class Target():

    def __init__(self, indir):
        self.filelist = glob.glob(f"{indir}*pdb")
        self.voxel_dict_all = {}
        self.freq_dict = {}

    def set_grid(self, center, radius):
        self.grid = gr.Grid(center, radius)
        self.grid.generate_voxels()

    def get_structure(self, file):
        structure = read_pdb(file)
        return structure

    def get_grid_atoms(self, structure):
        structure_atoms = list(structure.get_atoms())
        dist = np.sqrt(3)*self.grid.radius #get dist from center to vertex
        atoms_near = neighbor_search(structure_atoms, self.grid.center, dist)
        grid_atoms = [atom for atom in atoms_near \
                        if all(self.grid.v1 <= atom.get_coord()) \
                        and all(atom.get_coord() <= self.grid.v8)]
        return grid_atoms

    def check_voxels(self, grid_atoms):
        voxel_centers = np.array([v.center for v in self.grid.voxels])
        atom_coords = np.array([a.get_coord() for a in grid_atoms])
        dist = distance.cdist(atom_coords, voxel_centers, 'euclidean')
        min = dist.argmin(axis=1) #get index of closest voxel to each atom
        voxel_dict = {}
        for i, atom in enumerate(grid_atoms):
            voxel = min[i]
            element = atom.get_name()[0]
            voxel_dict.setdefault(voxel, []).append(element)
        return voxel_dict

    def analyze_trajectory(self, file):
        model = self.get_structure(file)
        grid_atoms = self.get_grid_atoms(model)
        voxel_dict = self.check_voxels(grid_atoms)
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

def read_pdb(file):
    parser = PDBParser()
    pdb_id, ext = file.split(".pdb")
    structure = parser.get_structure(pdb_id, file)
    return structure

def neighbor_search(atom_list, center, distance):
    neighbor_search = NeighborSearch(atom_list)
    atoms_near = neighbor_search.search(center, distance, 'A')
    return atoms_near


if __name__ == "__main__":
    target = Target("0/trajectory_1*")
    target.set_grid((2.173, 15.561, 28.257), 5)
    p = Pool(1)
    dicts = p.map(target.analyze_trajectory, target.filelist)
    p.close()
    p.join()
    for d in dicts:
        target.merge_voxel_dicts(d)
    print(target.voxel_dict_all)
    target.get_frequencies()


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
