import os
import numpy as np
import mdtraj as md
from multiprocessing import Pool
from functools import partial
from sklearn.neighbors import KDTree
# from Bio.PDB.PDBParser import PDBParser
# from Bio.PDB.NeighborSearch import NeighborSearch

def load_topology(file):
    return md.load(file).topology

def get_indices(topology, resname, atoms):
    if isinstance(atoms, list) or isinstance(atoms, tuple):
        atoms = " ".join(atoms)
    query = f"resname {resname} and name {atoms}"
    return topology.select(query)

def load_trajectory(file, indices=None):
    return md.load(file, atom_indices=indices)

def read_pdb(file):
    parser = PDBParser()
    pdb_id, ext = file.split(".pdb")
    structure = parser.get_structure(pdb_id, file)
    return structure

def neighbor_search_biopython(atom_list, center, distance):
    neighbor_search = NeighborSearch(atom_list)
    atoms_near = neighbor_search.search(center, distance, 'A')
    return atoms_near

def neighbor_search(coordinates, center, dist):
    center = np.array(center).reshape(1, 3)
    tree = KDTree(coordinates, leaf_size=3)
    result = tree.query_radius(center, r=dist)
    ind = result[0]
    atoms_near = coordinates[ind]
    return atoms_near

def format_line_pdb(coords, atomname, bfact, models = None, atomnum = "1", resname="UNK", chain="A", resnum="1", occ=1.00):
    x, y, z = coords
    atomstr = "ATOM"
    element = atomname[0]
    line = f"{atomstr:5}{atomnum:>5} {atomname:4} {resname:3} {chain}{resnum:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:6.2f}{bfact:7.2f}{element:>12}{models}\n"
    return line

def inside_grid(coord, lower_coord, upper_coord):
    return all(lower_coord <= coord) and all(coord <= upper_coord)

def basename_without_extension(filename):
    basename = os.path.basename(filename)
    basename, ext = os.path.splitext(basename)
    return basename

def frequency_dict(dict_, key, value):
    if dict_ is None:
        dict_ = {}
    if key in dict_:
        dict_[key] += value
    else:
        dict_[key] = value
    return dict_

def list_dict(dict_, key, value):
    if dict_ is None:
        dict_ = {}
    dict_.setdefault(key, []).append(value)
    return dict_


def merge_array_dicts(*dicts):
    merged_dict = {}
    union_keys = set().union(*dicts)

    for key in union_keys:
        lst = [i for i, d in enumerate(dicts) if key in d]
        merged_dict[key] = np.vstack([dicts[i][key] for i in lst])

    return merged_dict


def custom_path(dir, custom_var, string, ext):
    return os.path.join(dir, f"{custom_var}{string}{ext}")

def accepted_pele_steps(report):
    delimiter = " "*4
    with open(report) as r:
        columns = list(zip(*(line.strip().split(delimiter) for line in r)))

    for column in columns:
        if column[0] == 'numberOfAcceptedPeleSteps':
            accepted_steps = [int(value) for value in column[1:]]
    return accepted_steps


def parallelize(func, iterable, n_workers, **kwargs):
    f = partial(func, **kwargs)
    if n_workers > 1:
        with Pool(n_workers) as p:
            return p.map(f, iterable)
    else:
        return list(map(f, iterable))

def centroid(coords):
    return coords.sum(axis = 0) / len(coords)

def midpoint(point, other_point):
    x, y, z = point
    other_x, other_y, other_z = other_point
    mid = np.array(((x + other_x)/2 , (y + other_y)/2, (z + other_z)/2))
    return mid
