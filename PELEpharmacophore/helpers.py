import os
import numpy as np
import mdtraj as md
from multiprocessing import Pool
from functools import partial
from sklearn.neighbors import KDTree


def load_topology(file):
    return md.load(file).topology


def get_indices(topology, resname, atoms):
    if isinstance(atoms, list) or isinstance(atoms, tuple):
        atoms = " ".join(atoms)
    query = f"resname {resname} and name {atoms}"
    return topology.select(query)


def get_coordinates_from_trajectory(residue_name, trajectory, remove_hydrogen=False,
                                     only_first_model=False,
                                     indices_to_retrieve=None):
    """
    Given the path of a trajectory, it returns the array of coordinates
    that belong to the chosen residue.
    The resulting array will have as many elements as models the
    trajectory has (excluding the first model which will be always
    skipped).
    This method is prone to be parallelized.
    .. todo ::
       * Output warnings with logger, not with print.
    Parameters
    ----------
    residue_name : str
        The name of the residue whose coordinates will be extracted
    remove_hydrogen : bool
        Whether to remove all hydrogen atoms from the extracted
        coordinates array or not. Default is True
    trajectory : str
        The trajectory to extract the coordinates from
    only_first_model : bool
        Whether to retrieve the coordinates of the first model in the
        trajectory or all of them. It is optional and its default
        value is False
    indices_to_retrieve : list[int]
        The indices of the residue to extract. Default is None and
        will extract all of them
    Returns
    -------
    coordinates : a numpy.Array object
        The resulting array of coordinates
    """
    import numpy as np
    coordinates = list()
    # In case MODEL section is missing
    current_index = -1
    model_coords = []
    if only_first_model:
        skip_initial_structures = False
    else:
        skip_initial_structures = False
    with open(trajectory) as f:
        inside_model = False
        current_model = 0
        for i, line in enumerate(f):
            if len(line) <= 6:
                continue
            line_type = line[0:6]
            if line_type == "MODEL ":
                if inside_model:
                    print('Warning: ENDMDL declaration for model ' +
                          '{} might be missing'.format(current_model))
                inside_model = True
                current_model += 1
                current_index = -1
                model_coords = []
            if line_type == "ENDMDL":
                if not inside_model:
                    print('Warning: MODEL declaration for model ' +
                          '{} might be missing'.format(current_model + 1))
                inside_model = False
                # Only add the current model coordinates if the array is
                # not empty (to fulfill the dimensionality later on)
                if len(model_coords) > 0:
                    coordinates.append(np.array(model_coords))
                    model_coords = []
                # In case we are only interested in obtaining the
                # coordinates of the first model, we are done
                if only_first_model and current_model == 1:
                    break
            # First model will always be skipped, unless otherwise
            # established
            if current_model == 1 and skip_initial_structures:
                continue
            if line_type == "ATOM  " or line_type == "HETATM":
                current_residue_name = line[17:20]
                if current_residue_name == residue_name:
                    # Add one to current index (it initially equals -1)
                    current_index += 1
                    # In case we are interested in specific residue
                    # indices, retrieve only those
                    if (indices_to_retrieve is not None and
                            current_index not in indices_to_retrieve):
                        continue
                    # In case we have information about the element
                    # and we want to skip hydrogen atoms, do so
                    if remove_hydrogen and len(line) >= 78:
                        element = line[76:78]
                        element = element.strip()
                        element = element.strip(' ')
                        if element == 'H':
                            continue
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                    except ValueError:
                        print('Warning: invalid PDB format found in ' +
                              'line {}'.format(i) +
                              'of trajectory {}. '.format(trajectory) +
                              'Its coordinates will be skipped.')
                    point = np.array((x, y, z))
                    model_coords.append(point)
    # In case MODEL section was missing
    if not inside_model and len(model_coords) > 0:
        coordinates.append(np.array(model_coords))
    coordinates = np.array(coordinates)
    # When np.array.shape does not return a tuple of len 3 is because
    # its subarrays does not share the same dimensionality, so ligand
    # sizes are different.
    try:
        n_models_loaded, ligand_size, spatial_dimension = \
            coordinates.shape
    except ValueError:
        if len(coordinates) > 0:
            print('Warning: trajectory {} '.format(trajectory) +
                  'has an inconsistent ligand size throughout the ' +
                  'models. Its coordinates will be skipped.')
        # Return empty array
        return np.array(())
    if (n_models_loaded != current_model - 1 or spatial_dimension != 3) \
            and skip_initial_structures:
        print('Warning: unexpected dimensions found in the ' +
              'coordinate array from trajectory {}. '.format(trajectory) +
              'Its coordinates will be skipped.')
        # Return empty array
        return np.array(())

    return coordinates


def load_trajectory(file, indices=None):
    return md.load(file, atom_indices=indices)


def neighbor_search(coords, center, dist):
    center = np.array(center).reshape(1, 3)
    tree = KDTree(coords)
    result = tree.query_radius(center, r=dist)
    ind = result[0]
    atoms_near = coords[ind]
    return atoms_near


def format_line_pdb(coords, atomname, bfact=1, models = None, atomnum = "1", resname="UNK", chain="A", resnum="1", occ=1.00):
    x, y, z = coords
    atomstr = "HETATM"
    element = atomname[0]
    line = f"{atomstr:5}{atomnum:>5}{atomname:>4}{resname:>5}{chain:>2}{resnum:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:6.2f}{bfact:6.2f}{element:>12}\n"
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


def dict_function(func, dict_):
    return {key:func(value) for key, value in dict_}


def merge_array_dicts(*dicts):
    merged_dict = {}
    union_keys = set().union(*dicts)

    for key in union_keys:
        lst = [i for i, d in enumerate(dicts) if key in d] # ids of the dictionaries where certain key is present
        merged_dict[key] = np.vstack([dicts[i][key] for i in lst])

    return merged_dict

def gen_array_dicts(*dicts):
    gen_dict = {}
    union_keys = set().union(*dicts)

    for key in union_keys:
        gen_dict[key] = (d[key] for d in dicts if key in d)

    return gen_dict

def custom_path(dir, filename, ext):
    return os.path.join(dir, f"{filename}.{ext}")


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
            results = p.map(f, iterable)
    else:
        results = list(map(f, iterable))
    return results


def centroid(coords):
    return coords.sum(axis = 1) / coords.shape[1]


def pdbconvert(file, in_format="pdb", out_format="mae", outdir="."):
    input_path = os.path.abspath(file)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    convert_output = os.path.basename(file.replace(f".{in_format}", f".{out_format}"))
    convert_output = os.path.join(outdir, convert_output)
    schrodinger_path ="$SCHRODINGER/utilities/pdbconvert"
    command = f"{schrodinger_path} -i{in_format} {input_path} -o{out_format} {convert_output}"
    os.system(command)
