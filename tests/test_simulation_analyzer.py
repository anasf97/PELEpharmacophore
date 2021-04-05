import os
import numpy as np
import shutil
import pytest

DIR = os.path.dirname(__file__)
SIMULATION_1 = os.path.join(DIR, "data/simulation_1")
TRAJECTORY_1 = os.path.join(DIR, "data/simulation_1/output/0/trajectory_1.pdb")
TOPOLOGY = os.path.join(DIR, "data/simulation_1/output/topologies/topology_0.pdb")

EXPECTED_INDICES = {'HBD': np.array([5671]),
                    'HBA': np.array([5656, 5665, 5673]),
                    'ALI': np.array([5657, 5682]),
                    'ARO': np.array([5663, 5676])}

def test_get_indices(grid_analyzer, top_file=TOPOLOGY, expected_indices=EXPECTED_INDICES):
    ga = grid_analyzer
    topology = ga.get_topology(top_file)
    indices = ga.get_indices(topology, ga.resname)
    assert np.all(list(indices.values()) == list(expected_indices.values()))

EXPECTED_COORDS = [
                   [ 1.396, 14.473, 28.174],
                   [ 0.499,  9.416, 28.725],
                   [ 8.342, 11.331, 29.886],
                   [ 5.928, 12.535, 29.302],
                   [ 2.487, 14.577, 28.419],
                   [ 1.666, 22.194, 28.48 ],
                   [ 1.363, 18.291, 27.817],
                   [-0.923, 21.511, 27.99 ],
                   [ 3.414, 15.216, 28.589]
                  ]

def test_grid_atoms(grid_analyzer, atom_coords, trajectory=TRAJECTORY_1, expected_coords=EXPECTED_COORDS):
    sa = grid_analyzer
    topology = sa.get_topology(trajectory)

    grid_atoms = sa.get_grid_atoms(str[0])
    expected_coords = np.array([np.array(lst) for lst in expected_coords])

    assert atom_coords(grid_atoms).all() == expected_coords.all()


EXPECTED_VOXELS_CV = [77, 582, 736, 1070, 1281, 1295, 1854, 2400, 2639]

def test_check_voxels(grid_analyzer, active_voxels, trajectory=TRAJECTORY_1, expected_voxels=EXPECTED_VOXELS_CV):
    sa = grid_analyzer
    str = sa.get_structure(trajectory)

    grid_atoms = sa.get_grid_atoms(str[0])
    grid = sa.check_voxels(grid_atoms)

    assert active_voxels(grid) == expected_voxels



EXPECTED_VOXELS_MG = [77, 582, 736, 931, 1070, 1267, 1281, 1295, 1842, 1854, 2400, 2414, 2626, 2639]

def test_merge_grids(grid_analyzer, active_voxels, trajectory=TRAJECTORY_1, expected_voxels=EXPECTED_VOXELS_MG):
    sa = grid_analyzer
    str = sa.get_structure(trajectory)

    grid_atoms1 = sa.get_grid_atoms(str[0])
    grid_atoms2 = sa.get_grid_atoms(str[1])

    grid1 = sa.check_voxels(grid_atoms1)
    grid2 = sa.check_voxels(grid_atoms2)

    final_grid = sa.merge_grids(grid1, grid2)

    assert active_voxels(final_grid) == expected_voxels



EXPECTED_VOXELS_AT = [77, 582, 736, 931, 1070, 1267, 1281, 1295, 1842, 1854, 2400, 2414, 2626, 2639]

def test_analyze_trajectory_grid(grid_analyzer, active_voxels, expected_voxels=EXPECTED_VOXELS_AT):
    sa = grid_analyzer

    traj_grid = sa.analyze_trajectory(sa.traj_and_reports[0])

    assert active_voxels(traj_grid) == expected_voxels



EXPECTED_DICT = {'HBA': 1.1, 'ALI': 1.1, 'ARO': 1.1, 'HBD': 0.6}

def test_set_frequency_filter_grid(grid_analyzer, expected_dict=EXPECTED_DICT):
    sa = grid_analyzer
    sa.grid = sa.analyze_trajectory(sa.traj_and_reports[0])

    sa.set_frequency_filter(1)

    assert sa.threshold_dict == expected_dict



EXPECTED_PHARMACOPHORES_PATH = os.path.join(DIR, "data/pharmacophores_1")

def test_save_pharmacophores_grid(grid_analyzer, compare_files, expected_pharmacophores_path=EXPECTED_PHARMACOPHORES_PATH):
    sa = grid_analyzer
    sa.run(1)
    sa.set_frequency_filter(0)

    outdir = os.path.join(DIR, "result")
    sa.save_pharmacophores(outdir)

    result = [os.path.join(outdir, f) for f in os.listdir(outdir)]
    expected = [os.path.join(expected_pharmacophores_path, f) for f in os.listdir(expected_pharmacophores_path)]

    [compare_files(r, e) for r, e in zip(result, expected)]
    shutil.rmtree(outdir)



EXPECTED_COORDS_MS = [
                       [ 1.396, 14.473, 28.174],
                       [ 0.499,  9.416, 28.725],
                       [ 8.342, 11.331, 29.886],
                       [ 5.928, 12.535, 29.302],
                       [ 2.487, 14.577, 28.419],
                       [ 1.666, 22.194, 28.48 ],
                       [ 1.363, 18.291, 27.817],
                       [-0.923, 21.511, 27.99 ],
                       [ 3.414, 15.216, 28.589]
                     ]

def test_analyze_trajectory_meanshift(meanshift_analyzer, atom_coords, expected_coords=EXPECTED_COORDS_MS):
    sa = meanshift_analyzer

    featured_atoms = sa.analyze_trajectory(sa.traj_and_reports[0])
    expected_coords = np.array([np.array(lst) for lst in expected_coords])

    assert atom_coords(featured_atoms).all() == expected_coords.all()



EXPECTED_DICT_MS = {'HBD': 3.6, 'HBA': 2.2, 'ALI': 2.2, 'ARO': 2.2}

def test_set_frequency_filter_meanshift(meanshift_analyzer, active_voxels, expected_dict=EXPECTED_DICT_MS):
    sa = meanshift_analyzer
    sa.run(1)

    sa.set_frequency_filter(1)

    assert sa.threshold_dict == expected_dict


EXPECTED_PHARMACOPHORES_PATH_MS = os.path.join(DIR, "data/pharmacophores_ms")

def test_save_pharmacophores_meanshift(meanshift_analyzer, compare_files, expected_pharmacophores_path=EXPECTED_PHARMACOPHORES_PATH_MS):
    sa = meanshift_analyzer
    sa.run(1)
    sa.set_frequency_filter(0)

    outdir = os.path.join(DIR, "result_ms")
    sa.save_pharmacophores(outdir)

    result = [os.path.join(outdir, f) for f in os.listdir(outdir)]
    expected = [os.path.join(expected_pharmacophores_path, f) for f in os.listdir(expected_pharmacophores_path)]

    [compare_files(r, e) for r, e in zip(result, expected)]
    shutil.rmtree(outdir)
