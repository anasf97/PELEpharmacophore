import os
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch

def read_pdb(file):
    parser = PDBParser()
    pdb_id, ext = file.split(".pdb")
    structure = parser.get_structure(pdb_id, file)
    return structure

def neighbor_search(atom_list, center, distance):
    neighbor_search = NeighborSearch(atom_list)
    atoms_near = neighbor_search.search(center, distance, 'A')
    return atoms_near

def format_line_pdb(coords, atomname, bfact, models = None, atomnum = "1", resname="UNK", chain="A", resnum="1", occ=1.00):
    x, y, z = coords
    atomstr = "ATOM"
    element = atomname[0]
    line = f"{atomstr:5}{atomnum:>5} {atomname:4} {resname:3} {chain}{resnum:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:6.2f}{bfact:7.2f}{element:>12}{models}\n"
    return line

def atoms_inside_grid(atom, lower_coord, upper_coord):
    return all(lower_coord <= atom.get_coord()) and all(atom.get_coord() <= upper_coord)

def basename_without_extension(filename):
    basename = os.path.basename(filename)
    basename, ext = os.path.splitext(basename)
    return basename
