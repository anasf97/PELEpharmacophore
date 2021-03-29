import os
import numpy as np
from multiprocessing import Pool
from functools import partial
from sklearn.neighbors import KDTree
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch


def read_pdb(file):
    parser = PDBParser()
    pdb_id, ext = file.split(".pdb")
    structure = parser.get_structure(pdb_id, file)
    return structure

def neighbor_search_biopython(atom_list, center, distance):
    neighbor_search = NeighborSearch(atom_list)
    atoms_near = neighbor_search.search(center, distance, 'A')
    return atoms_near

def neighbor_search(atom_list, center, dist):
    coordinates = np.array([f.coordinates() for f in atom_list])
    center = np.array(center).reshape(1, 3)
    tree = KDTree(coordinates, leaf_size=3)
    ind = tree.query_radius(center, r=dist)
    ind = list(ind[0])
    atoms_near = [atom_list[i] for i in ind]
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
            output = pool.map(f, iterable)
    else:
        output = list(map(f, iterable))
        
    return output

def midpoint(point, other_point):
    x, y, z = point
    other_x, other_y, other_z = other_point
    mid = np.array(((x + other_x)/2 , (y + other_y)/2, (z + other_z)/2))
    return mid
